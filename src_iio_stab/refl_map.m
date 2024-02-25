% refl_map.m

clear all; close all;

PlotLambda = 1;
PlotRelative = 1;

RunNumber = 1;

if(RunNumber==1)
  N = 4; Nx = 16;
  Eps_Max = 1e-2; sigma = 1e-2;
elseif(RunNumber==2)
  N = 6; Nx = 24;
  Eps_Max = 0.1; sigma = 0.1;
elseif(RunNumber==3)
  N = 8; Nx = 32;
  Eps_Max = 0.1; sigma = 0.5;
elseif(RunNumber==100)
  % HOPS/AWE paper
  N = 16; Nx = 32;
  Eps_Max = 0.2; sigma = 0.99;
end
M=0;
Nz = 32;

alpha_bar = 0;
h_bar = 0.333; % What should this be?
gamma_u_bar = 1.21;
gamma_v_bar = 1.97;
gamma_w_bar = 2.23;
d = 2*pi;
c_0 = 1;
n_u = 1.0;
n_w = 1.1;   %2.3782, Carbon
% Mode = 1 (TE) or 2 (TM)
Mode = 2;
Taylor = true;

N_delta = 20;
N_Eps = 20;
%qq = [1,2];
qq = [1];
%qq = [1:6];

a = 1.0; b = 1.0;
% Nz = Nx
identy = eye(Nz+1);
[Dz,z] = cheb(Nz);

Eps = linspace(0,Eps_Max,N_Eps);

L = 2*pi;
k_u = sqrt(alpha_bar^2 + gamma_u_bar^2);
k_v = sqrt(alpha_bar^2 + gamma_v_bar^2);
k_w = sqrt(alpha_bar^2 + gamma_w_bar^2);
eta = (k_u+k_v+k_w)/3.0;
if(Mode==1)
  tau2 = 1.0;
  sigma2 = 1.0;
else
  tau2 = k_u^2/k_v^2;
  sigma2 = k_v^2/k_w^2;
end

% Create alphap, gamma_p, eep, and eem
[x,p,alphap,gamma_up,eep,eem] = setup_2d(Nx,L,alpha_bar,gamma_u_bar);
[x,p,alphap,gamma_vp,eep,eem] = setup_2d(Nx,L,alpha_bar,gamma_v_bar);
[x,p,alphap,gamma_wp,eep,eem] = setup_2d(Nx,L,alpha_bar,gamma_w_bar);

% Test functions
fu = cos(2*x);
fell = sin(2*x);
fu_x = real(ifft( (1i*p).*fft(fu) ));
fell_x = real(ifft( (1i*p).*fft(fell) ));


% Loop over q

R_taylor = zeros(N_Eps,N_delta);
R_pade = zeros(N_Eps,N_delta);
R_pade_safe = zeros(N_Eps,N_delta);

for s=1:length(qq)
  q = qq(s);
  if(N_delta==1)
    delta = [0];
  else
    delta = linspace(-sigma/(2*q+1),sigma/(2*q+1),N_delta);
  end
  
  omega_bar = c_0*(2*pi/d)*(q + 0.5);
  omega = (1+delta)*omega_bar;
  lambda = (2*pi*c_0./omega);

  tic;
  [U_n,Utilde_n,U,Utilde,V_u_n,Vtilde_u_n,V_u,Vtilde_u,...
      V_ell_n,Vtilde_ell_n,V_ell,Vtilde_ell,W_n,Wtilde_n,W,Wtilde] ...
      = mms_incidence2d(h_bar,eta,fu,fu_x,fell,fell_x,x,...
      alphap,gamma_up,gamma_vp,gamma_wp,Nx,N,Eps);

  if(Mode==1)
    zeta_n = (U_n - Utilde_n - V_u_n + Vtilde_u_n)/(-2*1i*eta);
    psi_n = (U_n + Utilde_n + V_u_n + Vtilde_u_n)/(-2);
    theta_n = (V_ell_n - Vtilde_ell_n - W_n + Wtilde_n)/(-2*1i*eta);
    mu_n = (V_ell_n + Vtilde_ell_n + W_n + Wtilde_n)/(-2);
  else
    zeta_n = (U_n - Utilde_n - V_u_n + Vtilde_u_n)/(-2*1i*eta);
    psi_n = (U_n + Utilde_n + tau2*V_u_n + tau2*Vtilde_u_n)/(-2);
    theta_n = (V_ell_n - Vtilde_ell_n - W_n + Wtilde_n)/(-2*1i*eta);
    mu_n = (V_ell_n + Vtilde_ell_n + sigma2*W_n + sigma2*Wtilde_n)/(-2);
  end

  [U_n_tfe,V_n_u_tfe,V_n_ell_tfe,W_n_tfe] ...
      = threelayer_iio_tfe_helmholtz(zeta_n,psi_n,theta_n,mu_n,...
      h_bar,eta,fu,fell,tau2,sigma2,p,alphap,gamma_up,gamma_vp,gamma_wp,...
      eep,eem,a,b,Nx,Nz,N);
  toc;

  % Correct order of dimensions - why do we need to do this?
  U_n_tfe = permute(U_n_tfe,[2 1]);
  V_n_u_tfe = permute(V_n_u_tfe,[2 1]);
  V_n_ell_tfe = permute(V_n_ell_tfe,[2 1]);
  W_n_tfe = permute(W_n_tfe,[2 1]);
  
  [ee_flat,ru_flat,rl_flat] ...
      = energy_defect(tau2,U_n_tfe,V_n_u_tfe,V_n_ell_tfe,W_n_tfe,...
      d,alpha_bar,gamma_u_bar,gamma_v_bar,gamma_w_bar,Eps,delta,...
      Nx,0,0,N_Eps,N_delta,1);
  
  % Don't compute pade_sum unless Taylor is false
  if Taylor == true
    [ee_taylor,ru_taylor,rl_taylor] ...
        = energy_defect(tau2,U_n_tfe,V_n_u_tfe,V_n_ell_tfe,W_n_tfe,...
        d,alpha_bar,gamma_u_bar,gamma_v_bar,gamma_w_bar,Eps,delta,...
        Nx,N,M,N_Eps,N_delta,1);
  else
    [ee_pade,ru_pade,rl_pade] ...
        = energy_defect(tau2,U_n_tfe,V_n_u_tfe,V_n_ell_tfe,W_n_tfe,...
        d,alpha_bar,gamma_u_bar,gamma_v_bar,gamma_w_bar,Eps,delta,...
        Nx,N,M,N_Eps,N_delta,2);
  end
  
  % We aren't currently testing pade_safe
  % [ee_pade_safe,ru_pade_safe,rl_pade_safe] ...
  %     = energy_defect(tau2,U_n_tfe,V_n_u_tfe,V_n_ell_tfe,W_n_tfe,...
  %     d,alpha_bar,gamma_u_bar,gamma_v_bar,gamma_w_bar,Eps,delta,...
  %     Nx,N,M,N_Eps,N_delta,3);
  
  % Plot the energy defect (log10)
  
  figure(1);
  set(gca,'FontSize',12);
  if(PlotLambda==0)
    if(Taylor == true)
      contourf(omega,Eps,log10(abs(ee_taylor))); hold on;
    else
      contourf(omega,Eps,log10(abs(ee_pade))); hold on;
    end
      xlabel('$\omega$','interpreter','latex','FontSize',18);
  else
    if(Taylor == true)
      contourf(lambda,Eps,log10(abs(ee_taylor))); hold on;
    else
      contourf(lambda,Eps,log10(abs(ee_pade))); hold on;
    end
    xlabel('$\lambda$','interpreter','latex','FontSize',18);
  end
  ylabel('$\varepsilon$','interpreter','latex','FontSize',20);
  title('$D$','interpreter','latex','FontSize',18);
  colorbar; colormap hot;


  figure(2);
  set(gca,'FontSize',12);
  if(PlotRelative==0)
    if(Taylor == true)
      RR = ru_taylor;
    else
      RR = ru_pade;
    end
  else
    if(Taylor == true)
      RR = ru_taylor./ru_flat;
    else
      RR = ru_pade./ru_flat;  
    end
  end
  if(PlotLambda==0)
    contourf(omega,Eps,RR); hold on;
    xlabel('$\omega$','interpreter','latex','FontSize',18);
  else
    contourf(lambda,Eps,RR); hold on;
    xlabel('$\lambda$','interpreter','latex','FontSize',18);
  end
  ylabel('$\varepsilon$','interpreter','latex','FontSize',20);
  title('$R$','interpreter','latex','FontSize',18);
  colorbar; colormap hot;
end
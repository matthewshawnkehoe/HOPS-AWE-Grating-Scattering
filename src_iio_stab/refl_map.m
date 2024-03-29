% refl_map.m

clear all; close all;

PlotLambda = 1;
PlotRelative = 1;

RunNumber = 100;

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
  % IIO Stability paper
  N = 16; Nx = 24;
  Eps_Max = 0.2; sigma = 0.99;
end
Ny = Nx;
Nz = 32;

alpha = 0;
beta = 0;
gamma_u = 1.21;
gamma_v = 1.97;
gamma_w = 2.23;
c_0 = 1;
n_u = 1.0;
n_w = 1.1;   %2.3782, Carbon
% Mode = 1 (TE) or 2 (TM)
Mode = 2;
Taylor = false;

N_Eps = 100;
%qq = [1,2];
qq = [1];
%qq = [1:6];

a = 0.5; b = 0.5;
identy = eye(Nz+1);
[Dz,z] = cheb(Nz);

Eps = linspace(0,Eps_Max,N_Eps);

d_x = 0.65;
d_y = 0.65;
k_u = sqrt(alpha^2 + gamma_u^2);
k_v = sqrt(alpha^2 + gamma_v^2);
k_w = sqrt(alpha^2 + gamma_w^2);
eta = (k_u+k_v+k_w)/3.0;
if(Mode==1)
  tau2 = 1.0;
  sigma2 = 1.0;
else
  tau2 = k_u^2/k_v^2;
  sigma2 = k_v^2/k_w^2;
end

% Create alphap, gamma_p, eep, and eem
[x,p,alphap,gamma_up,eep,eem] = setup_2d(Nx,d_x,alpha,gamma_u);
[x,p,alphap,gamma_vp,eep,eem] = setup_2d(Nx,d_x,alpha,gamma_v);
[x,p,alphap,gamma_wp,eep,eem] = setup_2d(Nx,d_x,alpha,gamma_w);
[y,p,alphap,gamma_up,eep,eem] = setup_2d(Ny,d_y,alpha,gamma_u);
[y,p,alphap,gamma_vp,eep,eem] = setup_2d(Ny,d_y,alpha,gamma_v);
[y,p,alphap,gamma_wp,eep,eem] = setup_2d(Ny,d_y,alpha,gamma_w);

% Test functions
fu = (1/4) * (cos(2*pi*x/d_x) + cos(2*pi*y/d_x));
fell = (1/4) * (cos(2*pi*x/d_x) + cos(2*pi*y/d_x));
fu_x = real(ifft( (1i*p).*fft(fu) ));
fell_x = real(ifft( (1i*p).*fft(fell) ));

% Loop over q
R_taylor = zeros(N_Eps,1);
R_pade = zeros(N_Eps,1);
R_pade_safe = zeros(N_Eps,1);

h_bar = pi/gamma_v + 10^(-16);

for s=1:length(qq)
  q = qq(s);
  omega = c_0*(2*pi/d_x)*(q + 0.5);
  % lambda = linspace(2*pi*c_0/omega+0.125,(2*pi*c_0./omega)+0.3,N_Eps)*1000;
  lambda = linspace(600, 750, N_Eps);

  k_u = n_u*omega/c_0;
  gamma_u_bar = sqrt(k_u^2 - alpha^2);
  [x,p,alphap,gamma_up,eep,eem] ...
      = setup_2d(Nx,d_x,alpha,gamma_u);

  k_w = n_w*omega/c_0;
  gamma_w = sqrt(k_w^2 - alpha^2);
  [x,p,alphap,gamma_wp,eep,eem] ...
      = setup_2d(Nx,d_x,alpha,gamma_w);

  tic;
  for N = 0:16
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
    
    % Correct order of dimensions - why do we need to do this?
    U_n_tfe = permute(U_n_tfe,[2 1]);
    V_n_u_tfe = permute(V_n_u_tfe,[2 1]);
    V_n_ell_tfe = permute(V_n_ell_tfe,[2 1]);
    W_n_tfe = permute(W_n_tfe,[2 1]);
      
    [ee_flat,ru_flat,rl_flat] ...
        = energy_defect(tau2,U_n_tfe,V_n_u_tfe,V_n_ell_tfe,W_n_tfe,...
        d_x,alpha,gamma_u,gamma_v,gamma_w,Eps,...
        Nx,0,N_Eps,1);
      
    % Don't compute pade_sum unless Taylor is false
    if Taylor == true
      [ee_taylor,ru_taylor,rl_taylor] ...
          = energy_defect(tau2,U_n_tfe,V_n_u_tfe,V_n_ell_tfe,W_n_tfe,...
          d_x,alpha,gamma_u,gamma_v,gamma_w,Eps,...
            Nx,N,N_Eps,1);
    else
      [ee_pade,ru_pade,rl_pade] ...
          = energy_defect(tau2,U_n_tfe,V_n_u_tfe,V_n_ell_tfe,W_n_tfe,...
          d_x,alpha,gamma_u,gamma_v,gamma_w,Eps,...
          Nx,N,N_Eps,2);
    end

    ee_pade_all(:,:,N+1) = ee_pade;
    ru_pade_all(:,N+1) = ru_pade;
    rl_pade_all(:,N+1) = rl_pade;
  end
  toc;
  % We aren't currently testing pade_safe
  % [ee_pade_safe,ru_pade_safe,rl_pade_safe] ...
  %     = energy_defect(tau2,U_n_tfe,V_n_u_tfe,V_n_ell_tfe,W_n_tfe,...
  %     d,alpha_bar,gamma_u_bar,gamma_v_bar,gamma_w_bar,Eps,...
  %     Nx,N,N_Eps,3);
  
  % Plot the energy defect (log10)
  
  % figure(1);
  % set(gca,'FontSize',12);
  % if(PlotLambda==0)
  %   if(Taylor == true)
  %     contourf(omega,Eps,log10(abs(ee_taylor))); hold on;
  %   else
  %     contourf(omega,Eps,log10(abs(ee_pade))); hold on;
  %   end
  %     xlabel('$\omega$','interpreter','latex','FontSize',18);
  % else
  %   if(Taylor == true)
  %     contourf(lambda,Eps,log10(abs(ee_taylor))); hold on;
  %   else
  %     contourf(lambda,Eps,log10(abs(ee_pade))); hold on;
  %   end
  %   xlabel('$\lambda$','interpreter','latex','FontSize',18);
  % end
  % ylabel('$\varepsilon$','interpreter','latex','FontSize',20);
  % title('$D$','interpreter','latex','FontSize',18);
  % colorbar; colormap hot;

  % 
  figure(2);
  set(gca,'FontSize',12);
  if(PlotRelative==0)
    if(Taylor == true)
      RR = ru_taylor_all;
    else
      RR = ru_pade_all;
    end
  else
    if(Taylor == true)
      RR = ru_taylor_all./ru_flat;
    else
      RR = ru_pade_all./ru_flat;  
    end
  end
  if(PlotLambda==0)
    contourf(omega,Eps,RR); hold on;
    xlabel('$\omega$ (nm)','interpreter','latex','FontSize',18);
  else
    lambda = repmat(lambda', 1, size(RR,2)); % How to fix this?
    Eps = repmat(Eps',1,size(RR,2));
    contourf(lambda,Eps,RR); hold on;
    xlabel('$\lambda$ (nm)','interpreter','latex','FontSize',18);
  end
  ylabel('$h$ (nm)','interpreter','latex','FontSize',20);
  title('$R$','interpreter','latex','FontSize',18);
  colorbar; colormap hot;
end
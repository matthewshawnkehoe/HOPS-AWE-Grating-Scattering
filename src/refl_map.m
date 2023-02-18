% refl_map.m

clear all; close all;

PlotLambda = 1;
PlotRelative = 1;

RunNumber = 100;

if(RunNumber==1)
  M = 4; Nx = 16;
  Eps_Max = 1e-2; sigma = 1e-2;
elseif(RunNumber==2)
  M = 6; Nx = 24;
  Eps_Max = 0.1; sigma = 0.1;
elseif(RunNumber==3)
  M = 8; Nx = 32;
  Eps_Max = 0.1; sigma = 0.5;
elseif(RunNumber==100)
  % HOPS/AWE paper
  M = 16; Nx = 32;
  Eps_Max = 0.2; sigma = 0.99;
end
N = M; Nz = 32;

alpha_bar = 0;
d = 2*pi;
c_0 = 1;
n_u = 1.0;
n_w = 1.1;   %2.3782, Carbon
% Mode = 1 (TE) or 2 (TM)
Mode = 2;
Taylor = true;

N_delta = 100;
N_Eps = 100;
%qq = [1,2];
%qq = [1];
qq = [1:6];

a = 1.0; b = 1.0;
% Nz = Nx
identy = eye(Nz+1);
[Dz,z] = cheb(Nz);

Eps = linspace(0,Eps_Max,N_Eps);

% Set up f, f_x

pp = (2*pi/d)*[0:Nx/2-1,-Nx/2:-1]';
xx = (d/Nx)*[0:Nx-1]';
%f =  (1/4)*sin(4*pi*xx/d);
%f_x = (pi/d)*cos(4*pi*xx/d);
%f = (1/2)*cos(2*pi*xx/d);
%f_x = -(2*pi/d)*sin(2*pi*xx/d);
%f = (1/4)*cos(4*xx);
%f_x = -sin(4*xx);
f = cos(xx);
f_x = -sin(xx);

% Sawtooth/Lipschitz
%P = 40;
%[f,f_x] = fourier_repr_lipschitz(P,xx);
% Rough Profile
%P = 40;
%[f,f_x] = fourier_repr_rough(P,xx);

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

  k_u_bar = n_u*omega_bar/c_0;
  gamma_u_bar = sqrt(k_u_bar^2 - alpha_bar^2);
  [xx,pp,alpha_bar_p,gamma_u_bar_p,eep,eem] ...
      = setup_2d(Nx,d,alpha_bar,gamma_u_bar);

  k_w_bar = n_w*omega_bar/c_0;
  gamma_w_bar = sqrt(k_w_bar^2 - alpha_bar^2);
  [xx,pp,alpha_bar_p,gamma_w_bar_p,eep,eem] ...
      = setup_2d(Nx,d,alpha_bar,gamma_w_bar);
  
  % MSK 09/21/22 Removed msk_
  [zeta_n_m,psi_n_m] = setup_zeta_psi_n_m(xx,pp,...
      alpha_bar,gamma_u_bar,f,f_x,Nx,N,M);
  
  if(Mode==1)
    tau2 = 1.0; % TE
  else
    tau2 = (n_u/n_w)^2; % TM
  end

  tic;
  [U_n_m,W_n_m,ubar_n_m,wbar_n_m] = ...
    two_layer_solve_fast(tau2,zeta_n_m,psi_n_m,...
    gamma_u_bar_p,gamma_w_bar_p,N,Nx,f,f_x,pp,alpha_bar,...
    gamma_u_bar,gamma_w_bar,Dz,a,b,Nz,M,identy,alpha_bar_p);
  toc;
  
  % Correct order of dimensions
  ubar_n_m_upd = permute(ubar_n_m,[3 2 1]);
  wbar_n_m_upd = permute(wbar_n_m,[3 2 1]);
  
  [ee_flat,ru_flat,rl_flat] ...
      = energy_defect(tau2,ubar_n_m_upd,wbar_n_m_upd,...
      d,alpha_bar,gamma_u_bar,gamma_w_bar,Eps,delta,...
      Nx,0,0,N_Eps,N_delta,1);
  
  % Don't compute pade_sum unless Taylor is false
  if Taylor == true
    [ee_taylor,ru_taylor,rl_taylor] ...
        = energy_defect(tau2,ubar_n_m_upd,wbar_n_m_upd,...
        d,alpha_bar,gamma_u_bar,gamma_w_bar,Eps,delta,...
        Nx,N,M,N_Eps,N_delta,1);
  else
    [ee_pade,ru_pade,rl_pade] ...
        = energy_defect(tau2,ubar_n_m_upd,wbar_n_m_upd,...
        d,alpha_bar,gamma_u_bar,gamma_w_bar,Eps,delta,...
        Nx,N,M,N_Eps,N_delta,2);
  end
  
  % We aren't currently testing pade_safe
  % [ee_pade_safe,ru_pade_safe,rl_pade_safe] ...
  %     = energy_defect(tau2,ubar_n_m_upd,wbar_n_m_upd,...
  %     d,alpha_bar,gamma_u_bar,gamma_w_bar,Eps,delta,...
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

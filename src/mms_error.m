% mms_error.m: Visualizes the error created by the method of manufactured
% solutions in the upper and lower fields.
%
% MSK 06/18/21 Created for testing when the lower layer has an upwards
% pointing normal derivative, +\partial_N w

% MSK 11/03/2021 - Updated to test section 6.5 of AWE paper

clear all; close all;

N = 20;
M = 20;
Nx = 1024;
N_delta = 100;
N_Eps = 100;
%qq = [0:5];
qq = [1];

%N = 8; M = 8;

%Eps_Max = 0.1;
Eps_Max = 2.0;
Eps = linspace(0,Eps_Max,N_Eps);
sigma = 0.99;
%sigma = 0.9;

alpha_bar = 0;
d = 2*pi;
c_0 = 1.0;
n_u = 1.0;
n_w = 1.1;
% Mode = 1 (TE) or 2 (TM)
Mode = 2;

A = 5.0;
B = 3.0;
r = 4;

a = 4.0;
b = 4.0;
Nz = 128;
identy = eye(Nz+1);
[Dz,z] = cheb(Nz);

% Set up f, f_x

pp = (2*pi/d)*[0:Nx/2-1,-Nx/2:-1]';
xx = (d/Nx)*[0:Nx-1]';
%f = exp(cos(xx));
%f_x = exp(cos(xx)).*(-sin(xx));
f = (1/4)*cos(4*xx);
f_x = -sin(4*xx);

% Lipschitz Profile
P = 120;
[f,f_x] = fourier_repr_lipschitz(P,xx);
% Rough Profile
% P = 120;
% [f,f_x] = fourier_repr_rough(P,xx);

% Loop over q

err_G_taylor = zeros(N_Eps,N_delta);
err_J_taylor = zeros(N_Eps,N_delta);
err_U_taylor = zeros(N_Eps,N_delta);
err_W_taylor = zeros(N_Eps,N_delta);
err_G_pade = zeros(N_Eps,N_delta);
err_J_pade = zeros(N_Eps,N_delta);
err_U_pade = zeros(N_Eps,N_delta);
err_W_pade = zeros(N_Eps,N_delta);
err_G_pade_safe = zeros(N_Eps,N_delta);
err_J_pade_safe = zeros(N_Eps,N_delta);
err_U_pade_safe = zeros(N_Eps,N_delta);
err_W_pade_safe = zeros(N_Eps,N_delta);
% BEGIN: DPN 6/25/21
err_ubar_taylor = zeros(N_Eps,N_delta);
err_wbar_taylor = zeros(N_Eps,N_delta);
err_ubar_pade = zeros(N_Eps,N_delta);
err_wbar_pade = zeros(N_Eps,N_delta);
err_ubar_pade_safe = zeros(N_Eps,N_delta);
err_wbar_pade_safe = zeros(N_Eps,N_delta);
% END: DPN 6/25/21

for s=1:length(qq)
  q = qq(s);
  delta = linspace(-sigma/(2*q+1),sigma/(2*q+1),N_delta);
  
  omega_bar = q+0.5;
  omega = (1+delta)*omega_bar;

  k_u_bar = n_u*omega_bar/c_0;
  gamma_u_bar = sqrt(k_u_bar^2 - alpha_bar^2);
  [xx,pp,alpha_bar_p,gamma_u_bar_p,eep,eem] ...
      = setup_2d(Nx,d,alpha_bar,gamma_u_bar);
  gamma_u_bar_r = gamma_u_bar_p(r+1);
  
  k_w_bar = n_w*omega_bar/c_0;
  gamma_w_bar = sqrt(k_w_bar^2 - alpha_bar^2);
  [xx,pp,alpha_bar_p,gamma_w_bar_p,eep,eem] ...
      = setup_2d(Nx,d,alpha_bar,gamma_w_bar);
  pp_r = pp(r+1);
  alpha_bar_r = alpha_bar_p(r+1);
  gamma_w_bar_r = gamma_w_bar_p(r+1);
  
  [xi_u_r_n_m,nu_u_r_n_m] = setup_xi_u_nu_u_n_m(A,r,xx,pp,...
      alpha_bar_p,gamma_u_bar_p,f,f_x,Nx,N,M);
  [u_n_m] = field_tfe_helmholtz_m_and_n(xi_u_r_n_m,f,pp,gamma_u_bar_p,...
      alpha_bar,gamma_u_bar,Dz,a,Nx,Nz,N,M,identy,alpha_bar_p);
  [G_n_m] = dno_tfe_helmholtz_m_and_n(u_n_m,f,pp,Dz,a,Nx,Nz,N,M);
  
  [xi_w_r_n_m,nu_w_r_n_m] = setup_xi_w_nu_w_n_m(B,r,xx,pp,...
      alpha_bar_p,gamma_w_bar_p,f,f_x,Nx,N,M);
  [w_n_m] = field_tfe_helmholtz_m_and_n_lf(xi_w_r_n_m,f,pp,gamma_w_bar_p,...
      alpha_bar,gamma_w_bar,Dz,b,Nx,Nz,N,M,identy,alpha_bar_p);
  [J_n_m] = dno_tfe_helmholtz_m_and_n_lf(w_n_m,f,pp,Dz,b,Nx,Nz,N,M);
  
  if(Mode==1)
    tau2 = 1.0; % TE
  else
    tau2 = (n_u/n_w)^2; % TM
  end

  zeta_r_n_m = xi_u_r_n_m - xi_w_r_n_m;
  psi_r_n_m = -nu_u_r_n_m - tau2*nu_w_r_n_m;
  
  [U_n_m,W_n_m,ubar_n_m,wbar_n_m] = ...
    two_layer_solve_fast(tau2,zeta_r_n_m,psi_r_n_m,gamma_u_bar_p,...
    gamma_w_bar_p,N,Nx,f,f_x,pp,alpha_bar,gamma_u_bar,...
    gamma_w_bar,Dz,a,b,Nz,M,identy,alpha_bar_p);
    
  nplot = zeros(N+1,1);
  mplot = zeros(M+1,1);
  
  for j=1:N_Eps
    for ell=1:N_delta
      alpha = (1+delta(ell))*alpha_bar;
      gamma_u = (1+delta(ell))*gamma_u_bar;
      gamma_w = (1+delta(ell))*gamma_w_bar;
      alpha_r = alpha_bar_r + delta(ell)*alpha_bar;
      k_u = (1+delta(ell))*k_u_bar;
      k_w = (1+delta(ell))*k_w_bar;
      gamma_u_r = sqrt(k_u^2 - alpha_r^2);
      gamma_w_r = sqrt(k_w^2 - alpha_r^2);
      
      % Exact
      xi_u_r = A*exp(1i*pp_r*xx).*exp(1i*gamma_u_r*Eps(j)*f);
      nu_u_r = (-1i*gamma_u_r + 1i*pp_r*Eps(j)*f_x).*xi_u_r;
      xi_w_r = B*exp(1i*pp_r*xx).*exp(-1i*gamma_w_r*Eps(j)*f);
      % Add a negative sign from taking -\partial_N w
      nu_w_r = (-1i*gamma_w_r - 1i*pp_r*Eps(j)*f_x).*xi_w_r;
      
      % BEGIN: DPN 6/25/21
      ubar = A*exp(1i*pp_r*xx).*exp(1i*gamma_u_r*a);
      wbar = B*exp(1i*pp_r*xx).*exp(-1i*gamma_w_r*(-b));
      % END: DPN 6/25/21

      % Taylor
      % G_approx = fcn_sum(1,G_n_m,Eps(j),delta(ell),Nx,N,M);
      % J_approx = fcn_sum(1,J_n_m,Eps(j),delta(ell),Nx,N,M);
      % U_approx = fcn_sum(1,U_n_m,Eps(j),delta(ell),Nx,N,M);
      % W_approx = fcn_sum(1,W_n_m,Eps(j),delta(ell),Nx,N,M);
      % % BEGIN: DPN 6/25/21
      % ubar_approx = fcn_sum(1,ubar_n_m,Eps(j),delta(ell),Nx,N,M);
      % wbar_approx = fcn_sum(1,wbar_n_m,Eps(j),delta(ell),Nx,N,M);
      % % END: DPN 6/25/21
      % err_G_taylor(j,ell) = norm(nu_u_r-G_approx,Inf);
      % err_J_taylor(j,ell) = norm(nu_w_r-J_approx,Inf);
      % err_U_taylor(j,ell) = norm(xi_u_r-U_approx,Inf);
      % err_W_taylor(j,ell) = norm(xi_w_r-W_approx,Inf);
      % % BEGIN: DPN 6/25/21
      % err_ubar_taylor(j,ell) = norm(ubar-ubar_approx,Inf);
      % err_wbar_taylor(j,ell) = norm(wbar-wbar_approx,Inf);
      % % END: DPN 6/25/21
      
      % Pade
      % G_approx = fcn_sum(2,G_n_m,Eps(j),delta(ell),Nx,N,M);
      % J_approx = fcn_sum(2,J_n_m,Eps(j),delta(ell),Nx,N,M);
      U_approx = fcn_sum(2,U_n_m,Eps(j),delta(ell),Nx,N,M);
      % W_approx = fcn_sum(2,W_n_m,Eps(j),delta(ell),Nx,N,M);
      ubar_approx = fcn_sum(2,ubar_n_m,Eps(j),delta(ell),Nx,N,M);
      % wbar_approx = fcn_sum(2,wbar_n_m,Eps(j),delta(ell),Nx,N,M);
      % err_G_pade(j,ell) = norm(nu_u_r-G_approx,Inf);
      % err_J_pade(j,ell) = norm(nu_w_r-J_approx,Inf);
      err_U_pade(j,ell) = norm(xi_u_r-U_approx,Inf);
      % err_W_pade(j,ell) = norm(xi_w_r-W_approx,Inf);
      err_ubar_pade(j,ell) = norm(ubar-ubar_approx,Inf);
      % err_wbar_pade(j,ell) = norm(wbar-wbar_approx,Inf);
      
      % Pade (safe)
      % G_approx = fcn_sum(3,G_n_m,Eps(j),delta(ell),Nx,N,M);
      % J_approx = fcn_sum(3,J_n_m,Eps(j),delta(ell),Nx,N,M);
      % U_approx = fcn_sum(3,U_n_m,Eps(j),delta(ell),Nx,N,M);
      % W_approx = fcn_sum(3,W_n_m,Eps(j),delta(ell),Nx,N,M);
      % BEGIN: DPN 6/25/21
      % ubar_approx = fcn_sum(3,ubar_n_m,Eps(j),delta(ell),Nx,N,M);
      % wbar_approx = fcn_sum(3,wbar_n_m,Eps(j),delta(ell),Nx,N,M);
      % END: DPN 6/25/21
      % err_G_pade_safe(j,ell) = norm(nu_u_r-G_approx,Inf);
      % err_J_pade_safe(j,ell) = norm(nu_w_r-J_approx,Inf);
      % err_U_pade_safe(j,ell) = norm(xi_u_r-U_approx,Inf);
      % err_W_pade_safe(j,ell) = norm(xi_w_r-W_approx,Inf);
      % BEGIN: DPN 6/25/21
      % err_ubar_pade_safe(j,ell) = norm(ubar-ubar_approx,Inf);
      % err_wbar_pade_safe(j,ell) = norm(wbar-wbar_approx,Inf);
      % END: DPN 6/25/21
    end
  end
  
  figure(1);

  contourf(omega,Eps,log10(err_U_pade)); hold on;
  set(gca,'FontSize',12)
  xlabel('$\omega=\underline{\omega}_1(1+\delta)$','interpreter','latex','FontSize',18);
  ylabel('$\varepsilon$','interpreter','latex','FontSize',20);
  title('Relative Error','interpreter','latex','FontSize',18);
  colorbar;

  figure(2);

  contourf(omega,Eps,log10(err_ubar_pade)); hold on;
  set(gca,'FontSize',12)
  xlabel('$\omega=\underline{\omega}_1(1+\delta)$','interpreter','latex','FontSize',18);
  ylabel('$\varepsilon$','interpreter','latex','FontSize',20);
  title('Relative Error','interpreter','latex','FontSize',18);
  colorbar;
end
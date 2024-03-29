% test_helmholtz_threelayer.m
%
% Script to test Helmholtz multiple layer solvers (2d).
%
% MSK 2/23/24

clear all; close all;

% Constants and values of Eps
Eps_values = [0.005, 0.01, 0.05, 0.1];
% Mode=1/2: TE/TM
Mode = 2;
L = 2*pi;
alpha = 0.1;
beta = 0.2;
hbar = 0.333;
gamma_u = 1.21;
gamma_v = 1.97;
gamma_w = 2.23;
Nx = 64;
Nz = 16;
a = 0.5;
b = 0.5;
N = 10;
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
[x,p,alphap,gamma_up,eep,eem] = setup_2d(Nx,L,alpha,beta);
[x,p,alphap,gamma_vp,eep,eem] = setup_2d(Nx,L,alpha,beta);
[x,p,alphap,gamma_wp,eep,eem] = setup_2d(Nx,L,alpha,beta);

% Test functions
fu = cos(2*x);
fell = sin(2*x);
fu_x = real(ifft( (1i*p).*fft(fu) ));
fell_x = real(ifft( (1i*p).*fft(fell) ));

fprintf('test_helmholtz_threelayer\n');

% Loop through all values of Eps and store data
for eps_idx = 1:length(Eps_values)
    Eps = Eps_values(eps_idx);
    fprintf('-------------------------\n');
    if(Mode==1)
      fprintf('TE Mode\n');
    else
      fprintf('TM Mode\n');
    end
    fprintf('alpha = %g  gamma_u = %g  gamma_v = %g  gamma_w = %g\n',...
        alpha,gamma_u,gamma_v,gamma_w);
    fprintf('hbar = %g\n',hbar);
    fprintf('Eps = %g  a = %g  b = %g\n',Eps,a,b);
    fprintf('Nx = %d  Nz = %d  N = %d\n',Nx,Nz,N);

    [U_n,Utilde_n,U,Utilde,V_u_n,Vtilde_u_n,V_u,Vtilde_u,...
        V_ell_n,Vtilde_ell_n,V_ell,Vtilde_ell,W_n,Wtilde_n,W,Wtilde] ...
        = mms_incidence2d(hbar,eta,fu,fu_x,fell,fell_x,x,...
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
        hbar,eta,fu,fell,tau2,sigma2,p,alphap,gamma_up,gamma_vp,gamma_wp,...
        eep,eem,a,b,Nx,Nz,N);
 
    [relerr_U,nplot_U] = compute_errors(U,U_n_tfe,Eps,N,Nx);

    filename = sprintf('three_%d_%g.mat',N,Eps);
    save(filename, 'U_n_tfe', 'V_n_u_tfe', 'V_n_ell_tfe', 'W_n_tfe', ....
         'relerr_U', 'nplot_U', 'Eps', 'N');
end
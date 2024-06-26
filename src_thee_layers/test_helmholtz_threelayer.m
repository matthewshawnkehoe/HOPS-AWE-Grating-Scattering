% test_helmholtz_threelayer.m
%
% Script to test Helmholtz multiple layer solvers (2d).
%
% MSK 2/23/24

% TODO - Adjust to your problem

clear all; close all;

% Constants and values of Eps
Eps_values = [0.005, 0.01, 0.05, 0.1];
N_delta = length(Eps_values);
% Mode=1/2: TE/TM
Mode = 2;
d = 2*pi;
alpha = 0.1;
%beta = 0.2;
hbar = 0.333;
gamma_u = 1.21;
gamma_v = 1.97;
gamma_w = 2.23;
Nx = 64;
Nz = 16;
a = 0.5;
b = 0.5;
N = 10;
M = 10;
k_u = sqrt(alpha^2 + gamma_u^2);
k_v = sqrt(alpha^2 + gamma_v^2);
k_w = sqrt(alpha^2 + gamma_w^2);
gamma_u = sqrt(k_u^2 - alpha^2);
gamma_v = sqrt(k_v^2 - alpha^2);
gamma_w = sqrt(k_w^2 - alpha^2);
eta = (k_u+k_v+k_w)/3.0;
if(Mode==1)
  tau2 = 1.0;
  sigma2 = 1.0;
else
  tau2 = k_u^2/k_v^2;
  sigma2 = k_v^2/k_w^2;
end

% Create alphap, gamma_p, eep, and eem
[x,p,alphap,gamma_up,eep,eem] = setup_2d(Nx,d,alpha,gamma_u);
[x,p,alphap,gamma_vp,eep,eem] = setup_2d(Nx,d,alpha,gamma_v);
[x,p,alphap,gamma_wp,eep,eem] = setup_2d(Nx,d,alpha,gamma_w);

% Test functions
fu = cos(2*x);
fell = sin(2*x);
fu_x = real(ifft( (1i*p).*fft(fu) ));
fell_x = real(ifft( (1i*p).*fft(fell) ));

% Sigma
sigma = sqrt(sigma2);

% Identity
identy = eye(Nz+1);

fprintf('test_helmholtz_threelayer\n');

% Loop through all values of Eps and Delta and store data
for delta_idx = 1:N_delta
    for eps_idx = 1:length(Eps_values)
        Eps = Eps_values(eps_idx);
        Delta_values  = linspace(-sigma,sigma,N_delta);
        Delta = Delta_values(delta_idx);

        fprintf('-------------------------\n');
        if(Mode==1)
          fprintf('TE Mode\n');
        else
          fprintf('TM Mode\n');
        end
        fprintf('alpha = %g  gamma_u = %g  gamma_v = %g  gamma_w = %g\n',...
            alpha,gamma_u,gamma_v,gamma_w);
        fprintf('hbar = %g\n',hbar);
        fprintf('Eps = %g  Delta = %g  a = %g  \n',Eps,Delta,a);
        fprintf('Nx = %d  Nz = %d  N = %d  M = %d\n',Nx,Nz,N,M);
    
        [U_nm,Utilde_nm,U,Utilde,V_u_nm,Vtilde_u_nm,V_u,Vtilde_u,...
            V_ell_nm,Vtilde_ell_nm,V_ell,Vtilde_ell,W_nm,Wtilde_nm,W,Wtilde] ...
            = mms_incidence2d(hbar,eta,fu,fu_x,fell,fell_x,x,...
            alphap,gamma_up,gamma_vp,gamma_wp,Nx,N,M,Eps,Delta);
    
        if(Mode==1)
          zeta_nm = (U_nm - Utilde_nm - V_u_nm + Vtilde_u_nm)/(-2*1i*eta);
          psi_nm = (U_nm + Utilde_nm + V_u_nm + Vtilde_u_nm)/(-2);
          theta_nm = (V_ell_nm - Vtilde_ell_nm - W_nm + Wtilde_nm)/(-2*1i*eta);
          mu_nm = (V_ell_nm + Vtilde_ell_nm + W_nm + Wtilde_nm)/(-2);
        else
          zeta_nm = (U_nm - Utilde_nm - V_u_nm + Vtilde_u_nm)/(-2*1i*eta);
          psi_nm = (U_nm + Utilde_nm + tau2*V_u_nm + tau2*Vtilde_u_nm)/(-2);
          theta_nm = (V_ell_nm - Vtilde_ell_nm - W_nm + Wtilde_nm)/(-2*1i*eta);
          mu_nm = (V_ell_nm + Vtilde_ell_nm + sigma2*W_nm + sigma2*Wtilde_nm)/(-2);
        end
    
        [U_nm_tfe,V_nm_u_tfe,V_nm_ell_tfe,W_nm_tfe] ...
            = threelayer_tfe_helmholtz(zeta_nm,psi_nm,theta_nm,mu_nm,...
            hbar,eta,fu,fell,tau2,sigma2,p,alpha,gamma_u,gamma_v,gamma_w,...
            alphap,gamma_up,gamma_vp,gamma_wp,eep,eem,a,b,Nx,Nz,N,M,identy);
     
        [relerr_U,nplot_U] = compute_errors(U,U_nm_tfe,Eps,N,Nx);
    
        % filename = sprintf('three_%d_%g.mat',N,Eps);
        % save(filename, 'U_n_tfe', 'V_n_u_tfe', 'V_n_ell_tfe', 'W_n_tfe', ....
        %      'relerr_U', 'nplot_U', 'Eps', 'N');
    end
end
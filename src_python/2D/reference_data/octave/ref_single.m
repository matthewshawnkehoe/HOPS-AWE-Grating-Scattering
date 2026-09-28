% Reference for test_single_eps_delta.m (RunNumber=1, DoTwoLayerTest=2), no plots
addpath('/home/claude/hops/src'); addpath('/home/claude/ref');
M = 8; Nx = 16; Eps = 1e-7; sigma = 1e-2; N = M+1; Nz = Nx;
q = 1; alpha_bar = 0; d = 2*pi; c_0 = 1.0; n_u = 1.0; n_w = 1.1; Mode = 2;
A_u = 3.0; A_w = 5.0; r = 2; a = 1.0; b = 1.0;
identy = eye(Nz+1); [Dz,z] = cheb(Nz);
xx = (d/Nx)*[0:Nx-1]'; f = cos(xx); f_x = -sin(xx);
delta = sigma/(2*q+1); omega_bar = q+0.5;
k_u_bar = n_u*omega_bar/c_0; gamma_u_bar = sqrt(k_u_bar^2 - alpha_bar^2);
[xx,pp,alpha_bar_p,gamma_u_bar_p,eep,eem] = setup_2d(Nx,d,alpha_bar,gamma_u_bar);
k_w_bar = n_w*omega_bar/c_0; gamma_w_bar = sqrt(k_w_bar^2 - alpha_bar^2);
[xx,pp,alpha_bar_p,gamma_w_bar_p,eep,eem] = setup_2d(Nx,d,alpha_bar,gamma_w_bar);
[xi_u_r_n_m,nu_u_r_n_m] = setup_xi_u_nu_u_n_m(A_u,r,xx,pp,alpha_bar_p,gamma_u_bar_p,f,f_x,Nx,N,M);
[u_n_m] = field_tfe_helmholtz_m_and_n(xi_u_r_n_m,f,pp,gamma_u_bar_p,alpha_bar,gamma_u_bar,Dz,a,Nx,Nz,N,M,identy,alpha_bar_p);
[G_n_m] = dno_tfe_helmholtz_m_and_n(u_n_m,f,pp,Dz,a,Nx,Nz,N,M);
[xi_w_r_n_m,nu_w_r_n_m] = setup_xi_w_nu_w_n_m(A_w,r,xx,pp,alpha_bar_p,gamma_w_bar_p,f,f_x,Nx,N,M);
[w_n_m] = field_tfe_helmholtz_m_and_n_lf(xi_w_r_n_m,f,pp,gamma_w_bar_p,alpha_bar,gamma_w_bar,Dz,b,Nx,Nz,N,M,identy,alpha_bar_p);
[J_n_m] = dno_tfe_helmholtz_m_and_n_lf(w_n_m,f,pp,Dz,b,Nx,Nz,N,M);
tau2 = (n_u/n_w)^2;
zeta_r_n_m = xi_u_r_n_m - xi_w_r_n_m; psi_r_n_m = -nu_u_r_n_m - tau2*nu_w_r_n_m;
[U_n_m,W_n_m,ubar_n_m,wbar_n_m] = two_layer_solve_fast(tau2,zeta_r_n_m,psi_r_n_m,gamma_u_bar_p,gamma_w_bar_p,N,Nx,f,f_x,pp,alpha_bar,gamma_u_bar,gamma_w_bar,Dz,a,b,Nz,M,identy,alpha_bar_p);
[U2,W2,ubar2,wbar2] = two_layer_solve(tau2,zeta_r_n_m,psi_r_n_m,gamma_u_bar_p,gamma_w_bar_p,N,Nx,f,f_x,pp,alpha_bar,gamma_u_bar,gamma_w_bar,Dz,a,b,Nz,M,identy,alpha_bar_p);
alpha_bar_r = alpha_bar_p(r+1); pp_r = pp(r+1);
alpha_r = alpha_bar_r + delta*alpha_bar; k_u = (1+delta)*k_u_bar; k_w = (1+delta)*k_w_bar;
gamma_u_r = sqrt(k_u^2 - alpha_r^2); gamma_w_r = sqrt(k_w^2 - alpha_r^2);
xi_u_r = A_u*exp(1i*pp_r*xx).*exp(1i*gamma_u_r*Eps*f);
xi_w_r = A_w*exp(1i*pp_r*xx).*exp(-1i*gamma_w_r*Eps*f);
nu_u_r = (-1i*gamma_u_r + 1i*pp_r*Eps*f_x).*xi_u_r;
nu_w_r = (-1i*gamma_w_r - 1i*pp_r*Eps*f_x).*xi_w_r;
ubar = A_u*exp(1i*pp_r*xx).*exp(1i*gamma_u_r*a);
err_U = zeros(N+1,M+1,3); err_G = err_U; err_ubar = err_U; err_u = err_U; err_W = err_U; err_J=err_U;
ll = [0:Nz]'; tilde_z = cos(pi*ll/Nz); z_prime = (a/2.0)*(tilde_z - 1.0) + a;
u = zeros(Nx,Nz+1);
for j=1:Nx, for ell=0:Nz, zz = (a-Eps*f(j))*z_prime(ell+1)/a + Eps*f(j); u(j,ell+1) = A_u*exp(1i*pp_r*xx(j)).*exp(1i*gamma_u_r*zz); end, end
for n=0:N, for m=0:M, for st=1:3
  err_U(n+1,m+1,st) = norm(xi_u_r-fcn_sum(st,U_n_m,Eps,delta,Nx,n,m),Inf);
  err_W(n+1,m+1,st) = norm(xi_w_r-fcn_sum(st,W_n_m,Eps,delta,Nx,n,m),Inf);
  err_G(n+1,m+1,st) = norm(nu_u_r-fcn_sum(st,G_n_m,Eps,delta,Nx,n,m),Inf);
  err_J(n+1,m+1,st) = norm(nu_w_r-fcn_sum(st,J_n_m,Eps,delta,Nx,n,m),Inf);
  err_ubar(n+1,m+1,st) = norm(ubar-fcn_sum(st,ubar_n_m,Eps,delta,Nx,n,m),Inf);
  err_u(n+1,m+1,st) = max(max(abs(u-vol_fcn_sum(st,u_n_m,Eps,delta,Nx,Nz,n,m))));
end, end, end
save('-v7','ref_single.mat','xi_u_r_n_m','nu_u_r_n_m','xi_w_r_n_m','nu_w_r_n_m','u_n_m','G_n_m','w_n_m','J_n_m','U_n_m','W_n_m','ubar_n_m','wbar_n_m','U2','W2','ubar2','wbar2','err_U','err_W','err_G','err_J','err_ubar','err_u');

% Reduced mms_error.m reference (same physics, fewer (eps,delta) samples)
addpath('/home/claude/hops/src'); addpath('/home/claude/ref');
N = 16; M = 16; Nx = 32; N_delta = 7; N_Eps = 6; q = 1;
Eps = linspace(0,0.2,N_Eps); sigma = 0.99;
alpha_bar = 0; d = 2*pi; c_0 = 1.0; n_u = 1.0; n_w = 1.1;
A = 5.0; B = 3.0; r = 4; a = 4.0; b = 4.0; Nz = 32;
identy = eye(Nz+1); [Dz,z] = cheb(Nz);
xx = (d/Nx)*[0:Nx-1]'; f = (1/4)*cos(4*xx); f_x = -sin(4*xx);
delta = linspace(-sigma/(2*q+1),sigma/(2*q+1),N_delta); omega_bar = q+0.5;
k_u_bar = n_u*omega_bar/c_0; gamma_u_bar = sqrt(k_u_bar^2 - alpha_bar^2);
[xx,pp,alpha_bar_p,gamma_u_bar_p,eep,eem] = setup_2d(Nx,d,alpha_bar,gamma_u_bar);
k_w_bar = n_w*omega_bar/c_0; gamma_w_bar = sqrt(k_w_bar^2 - alpha_bar^2);
[xx,pp,alpha_bar_p,gamma_w_bar_p,eep,eem] = setup_2d(Nx,d,alpha_bar,gamma_w_bar);
pp_r = pp(r+1); alpha_bar_r = alpha_bar_p(r+1);
[xi_u_r_n_m,nu_u_r_n_m] = setup_xi_u_nu_u_n_m(A,r,xx,pp,alpha_bar_p,gamma_u_bar_p,f,f_x,Nx,N,M);
[xi_w_r_n_m,nu_w_r_n_m] = setup_xi_w_nu_w_n_m(B,r,xx,pp,alpha_bar_p,gamma_w_bar_p,f,f_x,Nx,N,M);
tau2 = (n_u/n_w)^2;
zeta_r_n_m = xi_u_r_n_m - xi_w_r_n_m; psi_r_n_m = -nu_u_r_n_m - tau2*nu_w_r_n_m;
tic;
[U_n_m,W_n_m,ubar_n_m,wbar_n_m] = two_layer_solve_fast(tau2,zeta_r_n_m,psi_r_n_m,gamma_u_bar_p,gamma_w_bar_p,N,Nx,f,f_x,pp,alpha_bar,gamma_u_bar,gamma_w_bar,Dz,a,b,Nz,M,identy,alpha_bar_p);
toc
errU = zeros(N_Eps,N_delta,3); errubar = errU; errW = errU;
for j=1:N_Eps, for ell=1:N_delta
  alpha_r = alpha_bar_r + delta(ell)*alpha_bar; k_u = (1+delta(ell))*k_u_bar; k_w = (1+delta(ell))*k_w_bar;
  gamma_u_r = sqrt(k_u^2 - alpha_r^2); gamma_w_r = sqrt(k_w^2 - alpha_r^2);
  xi_u_r = A*exp(1i*pp_r*xx).*exp(1i*gamma_u_r*Eps(j)*f);
  xi_w_r = B*exp(1i*pp_r*xx).*exp(-1i*gamma_w_r*Eps(j)*f);
  ubar = A*exp(1i*pp_r*xx).*exp(1i*gamma_u_r*a);
  for st=1:3
    errU(j,ell,st) = norm(xi_u_r-fcn_sum(st,U_n_m,Eps(j),delta(ell),Nx,N,M),Inf);
    errW(j,ell,st) = norm(xi_w_r-fcn_sum(st,W_n_m,Eps(j),delta(ell),Nx,N,M),Inf);
    errubar(j,ell,st) = norm(ubar-fcn_sum(st,ubar_n_m,Eps(j),delta(ell),Nx,N,M),Inf);
  end
end, end
save('-v7','ref_mms.mat','U_n_m','W_n_m','ubar_n_m','wbar_n_m','errU','errW','errubar','Eps','delta');

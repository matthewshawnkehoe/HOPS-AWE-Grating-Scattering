% Full paper Fig. 9 / plots/test_scenarios.m run (dielectric, N=M=16, Taylor, 6 bands, 100 x 100)
% with the ORIGINAL MATLAB sources. ~25 min per band in Octave (energy_defect loop dominates).
addpath('/home/claude/hops/src'); addpath('/home/claude/ref');
M = 16; Nx = 32; Eps_Max = 0.2; sigma = 0.99; N = M; Nz = 32;
alpha_bar = 0; d = 2*pi; c_0 = 1; n_u = 1.0; n_w = 1.1;
N_delta = 100; N_Eps = 100; qq = [1:6];
a = 1.0; b = 1.0; identy = eye(Nz+1); [Dz,z] = cheb(Nz);
Eps = linspace(0,Eps_Max,N_Eps);
xx = (d/Nx)*[0:Nx-1]'; f = cos(xx); f_x = -sin(xx);
out = struct();
for s=1:length(qq)
  q = qq(s);
  delta = linspace(-sigma/(2*q+1),sigma/(2*q+1),N_delta);
  omega_bar = c_0*(2*pi/d)*(q + 0.5);
  k_u_bar = n_u*omega_bar/c_0; gamma_u_bar = sqrt(k_u_bar^2 - alpha_bar^2);
  [xx,pp,alpha_bar_p,gamma_u_bar_p,eep,eem] = setup_2d(Nx,d,alpha_bar,gamma_u_bar);
  k_w_bar = n_w*omega_bar/c_0; gamma_w_bar = sqrt(k_w_bar^2 - alpha_bar^2);
  [xx,pp,alpha_bar_p,gamma_w_bar_p,eep,eem] = setup_2d(Nx,d,alpha_bar,gamma_w_bar);
  [zeta_n_m,psi_n_m] = setup_zeta_psi_n_m(xx,pp,alpha_bar,gamma_u_bar,f,f_x,Nx,N,M);
  tau2 = (n_u/n_w)^2;
  tic;
  [U_n_m,W_n_m,ubar_n_m,wbar_n_m] = two_layer_solve_fast(tau2,zeta_n_m,psi_n_m,gamma_u_bar_p,gamma_w_bar_p,N,Nx,f,f_x,pp,alpha_bar,gamma_u_bar,gamma_w_bar,Dz,a,b,Nz,M,identy,alpha_bar_p);
  t1 = toc; tic;
  ub = permute(ubar_n_m,[3 2 1]); wb = permute(wbar_n_m,[3 2 1]);
  [ee_flat,ru_flat,rl_flat] = energy_defect(tau2,ub,wb,d,alpha_bar,gamma_u_bar,gamma_w_bar,Eps,delta,Nx,0,0,N_Eps,N_delta,2);
  [ee_t,ru_t,rl_t] = energy_defect(tau2,ub,wb,d,alpha_bar,gamma_u_bar,gamma_w_bar,Eps,delta,Nx,N,M,N_Eps,N_delta,1);
  printf('q=%d solve %.0f s, energy %.0f s\n', q, t1, toc); fflush(stdout);
  out.(sprintf('q%d',q)) = struct('ubar',ubar_n_m,'wbar',wbar_n_m,'ee_t',ee_t,'ru_t',ru_t,'rl_t',rl_t,'ru_flat',ru_flat,'ee_flat',ee_flat,'delta',delta,'Eps',Eps);
  save('-v7',sprintf('ref_refl_dielectric_full.mat'),'out');
end

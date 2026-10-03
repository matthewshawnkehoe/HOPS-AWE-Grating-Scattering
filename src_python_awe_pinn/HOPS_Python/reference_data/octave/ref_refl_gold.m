% Reduced refl_map.m references: silver (refl_map.m) and dielectric (test_scenarios.m)
addpath('/home/claude/hops/src'); addpath('/home/claude/ref');
cases = {'gold'};
for cc=1:1
M=15; n_w = 1.48 + 1.883*1i; fcase=4;
Nx = 32; Eps_Max = 0.2; sigma = 0.99; N = M; Nz = 32;
alpha_bar = 0; d = 2*pi; c_0 = 1; n_u = 1.0;
N_delta = 5; N_Eps = 4; qq = [1 3];
a = 1.0; b = 1.0; identy = eye(Nz+1); [Dz,z] = cheb(Nz);
Eps = linspace(0,Eps_Max,N_Eps);
xx = (d/Nx)*[0:Nx-1]';
if fcase==4, f = cos(4*xx); f_x = -4*sin(4*xx); else, f = cos(xx); f_x = -sin(xx); end
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
  toc
  ub = permute(ubar_n_m,[3 2 1]); wb = permute(wbar_n_m,[3 2 1]);
  [ee_flat,ru_flat,rl_flat] = energy_defect(tau2,ub,wb,d,alpha_bar,gamma_u_bar,gamma_w_bar,Eps,delta,Nx,0,0,N_Eps,N_delta,1);
  [ee_t,ru_t,rl_t] = energy_defect(tau2,ub,wb,d,alpha_bar,gamma_u_bar,gamma_w_bar,Eps,delta,Nx,N,M,N_Eps,N_delta,1);
  [ee_p,ru_p,rl_p] = energy_defect(tau2,ub,wb,d,alpha_bar,gamma_u_bar,gamma_w_bar,Eps,delta,Nx,N,M,N_Eps,N_delta,2);
  [ee_s,ru_s,rl_s] = energy_defect(tau2,ub,wb,d,alpha_bar,gamma_u_bar,gamma_w_bar,Eps,delta,Nx,N,M,N_Eps,N_delta,3);
  out.(sprintf('q%d',q)) = struct('zeta',zeta_n_m,'psi',psi_n_m,'U',U_n_m,'W',W_n_m,'ubar',ubar_n_m,'wbar',wbar_n_m, ...
     'ee_flat',ee_flat,'ru_flat',ru_flat,'ee_t',ee_t,'ru_t',ru_t,'rl_t',rl_t,'ee_p',ee_p,'ru_p',ru_p,'rl_p',rl_p,'ee_s',ee_s,'ru_s',ru_s);
end
save('-v7',sprintf('ref_refl_%s.mat',cases{cc}),'out');
end

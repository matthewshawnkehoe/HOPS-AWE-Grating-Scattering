function [ee,ru,rl] = energy_defect(tau2,ubar_n,vbar_n_u, vbar_n_ell,...
    wbar_n, d,alpha_bar,gamma_u_bar,gamma_v_bar,gamma_w_bar,Eps,...
    Nx,N,N_Eps,SumType)

k_u = sqrt(alpha_bar^2 + gamma_u_bar^2);
k_v = sqrt(alpha_bar^2 + gamma_v_bar^2);
k_w = sqrt(alpha_bar^2 + gamma_w_bar^2);
ru = zeros(N_Eps,1);
rl = zeros(N_Eps,1);
PropMode_u = zeros(Nx,1);
PropMode_v = zeros(Nx,1);
PropMode_w = zeros(Nx,1);


for j=1:N_Eps
  [xx,pp,alpha_p,gamma_u_p,eep,eem] = setup_2d(Nx,d,alpha_bar,gamma_u_bar);
  [xx,pp,alpha_p,gamma_v_p,eep,eem] = setup_2d(Nx,d,alpha_bar,gamma_v_bar);
  [xx,pp,alpha_p,gamma_w_p,eep,eem] = setup_2d(Nx,d,alpha_bar,gamma_w_bar);
    
  gamma_i = gamma_u_p(0+1);
    
  [ubar] = fcn_sum_single(SumType,ubar_n,Eps(j),Nx,N);
  [vbaru] = fcn_sum_single(SumType,vbar_n_u,Eps(j),Nx,N);
  [vbarell] = fcn_sum_single(SumType,vbar_n_ell,Eps(j),Nx,N);
  [wbar] = fcn_sum_single(SumType,wbar_n,Eps(j),Nx,N);
    
  ubarhat = fft(ubar)/Nx;
  vbaruhat = fft(vbaru)/Nx;
  vbarellhat = fft(vbarell)/Nx;
  wbarhat = fft(wbar)/Nx;
    
  alpha_p_2 = alpha_p.^2;
  k_u_2 = k_u^2;
  k_v_2 = k_v^2;
  k_w_2 = k_w^2;
    
  B = gamma_u_p./gamma_i.*abs(ubarhat).^2;
  C = gamma_v_p./gamma_i.*abs(vbaruhat).^2;
  D = gamma_v_p./gamma_i.*abs(vbarellhat).^2;
  E = gamma_w_p./gamma_i.*abs(wbarhat).^2;
    
  PropMode_u(alpha_p_2 < k_u_2) = 1;
  PropMode_v(alpha_p_2 < k_v_2) = 1;
  PropMode_w(alpha_p_2 < k_w_2) = 1;
    
  ru(j) = sum(PropMode_u.*B);
  rv(j) = sum(PropMode_v.*C) + sum(PropMode_v.*D);  % What should this be?
  rl(j) = tau2*sum(PropMode_w.*D);
end


ee = 1.0 - ru - rl;
%ee = 1.0 - ru - rv - rl;

return;
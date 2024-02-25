function [ee,ru,rl] = energy_defect(tau2,ubar_n_m,vbar_n_u_m, vbar_n_ell_m,...
    wbar_n_m, d,alpha_bar,gamma_u_bar,gamma_v_bar,gamma_w_bar,Eps,delta,...
    Nx,N,M,N_Eps,N_delta,SumType)

ru = zeros(N_Eps,N_delta);
rl = zeros(N_Eps,N_delta);

k_u_bar = sqrt(alpha_bar^2 + gamma_u_bar^2);
k_v_bar = sqrt(alpha_bar^2 + gamma_v_bar^2);
k_w_bar = sqrt(alpha_bar^2 + gamma_w_bar^2);

PropMode_u = zeros(Nx,1);
PropMode_v = zeros(Nx,1);
PropMode_w = zeros(Nx,1);

for ell=1:N_delta
  for j=1:N_Eps
    alpha = (1+delta(ell))*alpha_bar;
    gamma_u = (1+delta(ell))*gamma_u_bar;
    gamma_v = (1+delta(ell))*gamma_v_bar;
    gamma_w = (1+delta(ell))*gamma_w_bar;
    k_u = (1+delta(ell))*k_u_bar;
    k_v = (1+delta(ell))*k_v_bar;
    k_w = (1+delta(ell))*k_w_bar;
  
    [xx,pp,alpha_p,gamma_u_p,eep,eem] = setup_2d(Nx,d,alpha,gamma_u);
    [xx,pp,alpha_p,gamma_v_p,eep,eem] = setup_2d(Nx,d,alpha,gamma_v);
    [xx,pp,alpha_p,gamma_w_p,eep,eem] = setup_2d(Nx,d,alpha,gamma_w);
  
    gamma_i = gamma_u_p(0+1);
    
    [ubar] = fcn_sum_fast(SumType,ubar_n_m,Eps(j),delta(ell),Nx,N,M);
    [vbaru] = fcn_sum_fast(SumType,vbar_n_u_m,Eps(j),delta(ell),Nx,N,M);
    [vbarell] = fcn_sum_fast(SumType,vbar_n_ell_m,Eps(j),delta(ell),Nx,N,M);
    [wbar] = fcn_sum_fast(SumType,wbar_n_m,Eps(j),delta(ell),Nx,N,M);
    
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
    
    ru(j,ell) = sum(PropMode_u.*B);
    rv(j,ell) = sum(PropMode_v.*C) + sum(PropMode_v.*D);  % What should this be?
    rl(j,ell) = tau2*sum(PropMode_w.*D);
  end
end

%ee = 1.0 - ru - rl;
ee = 1.0 - ru - rv - rl;

return;
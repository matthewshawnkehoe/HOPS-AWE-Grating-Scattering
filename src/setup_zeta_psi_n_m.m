function [zeta_n_m,psi_n_m] = setup_zeta_psi_n_m(xx,pp,...
    alpha_bar,gamma_u_bar,f,f_x,Nx,N,M)

% MSK 7/26/21: Changed the size of zeta_n_m,psi_n_m from (Nx,N+1,M+1) to (Nx,M+1,N+1)

gamma_u_m = zeros(1,M+1);
gamma_u_m(0+1) = gamma_u_bar;
gamma_u_m(1+1) = gamma_u_bar;

alpha_u_m = zeros(1,M+1);
alpha_u_m(0+1) = alpha_bar;
alpha_u_m(1+1) = alpha_bar;

E_n_m = E_exp_lf(gamma_u_m,f,N,M);

zeta_n_m = zeros(Nx,M+1,N+1);
psi_n_m = zeros(Nx,M+1,N+1);

for n=0:N
  for m=0:M
    zeta_n_m(:,m+1,n+1) = -E_n_m(:,m+1,n+1);
    for ell=0:m
      psi_n_m(:,m+1,n+1) = psi_n_m(:,m+1,n+1) ...
          + (1i*gamma_u_m(m-ell+1)).*E_n_m(:,ell+1,n+1);
    end
    if(n>0)
      psi_n_m(:,m+1,n+1) = psi_n_m(:,m+1,n+1)...
          + f_x.*(1i*alpha_u_m(m-ell+1)).*E_n_m(:,ell+1,n-1+1);
    end    
  end
end

return;
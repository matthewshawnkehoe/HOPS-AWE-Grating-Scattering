function [E_nm_lf] = E_exp_lf(gamma_m,f,N,M)

Nx = length(f);
E_nm_lf = zeros(Nx,M+1,N+1);

% n = 0, m = 0
E_nm_lf(:,0+1,0+1) = ones(Nx,1);

% n = 0, m > 0: Do nothing as E_nm_lf(:,0+1,m+1) = 0

% n > 0, m >= 0

for n=1:N
  for m=0:M
    sum = zeros(Nx,1);
    for r=0:m
      sum = sum + E_nm_lf(:,r+1,n-1+1)*(1i*gamma_m(m-r+1));
    end
    E_nm_lf(:,m+1,n+1) = -f.*sum/n;
  end
end

return;
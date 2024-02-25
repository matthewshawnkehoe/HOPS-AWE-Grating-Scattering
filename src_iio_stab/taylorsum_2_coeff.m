function [tsum] = taylorsum_2_coeff(c,Eps,delta,N,M)
% taylorsum - Sums a truncated Taylor series.
%
% Inputs:
%
% c - Taylor series coefficients: [c_0,...,c_N]
% Eps - Value at which to sum for N
% Delta - Value at which to sum for M
% N - Degree of truncated Taylor series
% M - Degree of truncated Taylor series
%
% Outputs:
%
% tsum - Taylor sum evaluated at Eps and Delta
%
% DPN 2/7/12
% MSK 7/21/21 Optimized double for loop

rho = sqrt(Eps^2 + delta^2);
theta = atan2(delta,Eps);

N_M_min = min(N,M);
coeff = zeros(N_M_min+1,1);
c1 = cos(theta).^(0:N_M_min);
c2 = sin(theta).^(0:N_M_min);
for p=0:N_M_min
  for q=0:p
    %ii = index_nm(p-q,q,N,M);
    coeff(p+1) = coeff(p+1) + c(p-q+1,q+1)*c1(p-q+1)*c2(q+1);
  end
end

N_M_min_over_2 = floor(N_M_min/2.0);
tsum = taylorsum(coeff,rho,N_M_min_over_2);

return;
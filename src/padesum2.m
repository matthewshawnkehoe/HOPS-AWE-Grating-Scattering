function [psum] = padesum2(c,Eps,delta,N,M)
% padesum2 - Uses Pade approximation to sum a truncated Taylor series.
%
% Inputs:
%   c - Taylor series coefficients: [c_{0,0},...,c_{0,N},...,c_{M,0},...,c_{M,N}]
%   Eps, delta - Value at which to sum
%   N,M - Degree of truncated Taylor series
%
% Output:
%   psum: Taylor sum evaluated at Eps, delta
%
% DPN 12/20/15
% MSK 7/21/21 Optimized double for loop

rho = sqrt(Eps^2 + delta^2);
theta = atan2(delta,Eps);

% Form the ctilde

N_M_min = min(N,M);

coeff = zeros(N_M_min+1,1);
c1 = cos(theta).^(0:N_M_min);
c2 = sin(theta).^(0:N_M_min);
for p=0:N_M_min
  for q=0:p
    coeff(p+1) = coeff(p+1) + c(p-q+1,q+1)*c1(p-q+1)*c2(q+1);
  end
end

N_M_min_over_2 = floor(N_M_min/2.0);
psum = padesum(coeff,rho,N_M_min_over_2);

return;
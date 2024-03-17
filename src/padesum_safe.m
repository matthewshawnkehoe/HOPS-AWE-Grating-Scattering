function [psum,a,b] = padesum_safe(c,Eps,M)
% padesum_safe.m: Sums a truncated Taylor series via Pade approximation safely.
%
% Inputs:
%   c: Taylor series coefficients: [c_0,...,c_N]
%   Eps: Value at which to sum
%   M: Numerator and denominator degree (M = N/2)
%
% Outputs:
%   psum: Pade approximant evaluated at Eps
%   a: Numerator coefficients: [a_0,...,a_M]
%   b: Denominator coefficients: [1,b_1,...,b_M]
%
% Note: Assumes that N = 2*M (otherwise ignores c_N)
%
% DPN 2/7/12
% DPN 11/25/14 (safe version)

N_true = -1;
for j=1:2*M+1
  if(abs(c(j))>1e-14)
    N_true = j;
  end
end
M_safe = floor(N_true/2.0);

if(M_safe==0)
  a = c(1);
  b = [1];
else
  H = toeplitz(c(M_safe+1:2*M_safe-1+1),c(M_safe+1:-1:1+1));
  bb = -H\c(M_safe+1+1:2*M_safe+1);
  b = [1;bb];
  aa = conv(b,c(0+1:M_safe+1));
  a = aa(0+1:M_safe+1);
end
psum = polyval(a(end:-1:1),Eps)./polyval(b(end:-1:1),Eps);
function [psum,a,b] = padesum(c,Eps,M)
% padesum - Sums a truncated Taylor series via Pade approximation.
%
% Inputs:
%
% c - Taylor series coefficients: [c_0,...,c_N]
% Eps - Value at which to sum
% M - Numerator and denominator degree (M = N/2)
%
% Outputs:
%
% psum - Pade approximant evaluated at Eps
% a - Numerator coefficients: [a_0,...,a_M]
% b - Denominator coefficients: [1,b_1,...,b_M]
%
% Note: Assumes that N = 2*M (otherwise ignores c_N)
%
% DPN 2/7/12

if(M==0)
  a = c(1);
  b = [1];
else
  H = toeplitz(c(M+1:2*M-1+1),c(M+1:-1:1+1));
  bb = -H\c(M+1+1:2*M+1);
  b = [1;bb];
  aa = conv(b,c(0+1:M+1));
  a = aa(0+1:M+1);
end
psum = polyval(a(end:-1:1),Eps)./polyval(b(end:-1:1),Eps);
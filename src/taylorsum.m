function [tsum] = taylorsum(c,Eps,N)
% taylorsum.m: Sums a truncated Taylor series.
%
% Inputs:
%   c: Taylor series coefficients for [c_0,...,c_N]
%   Eps: Value at which to sum
%   N: Degree of truncated Taylor series
%
% Output:
%   tsum: Taylor sum evaluated at Eps
%
% DPN 2/7/12
% MSK 2/21/21 - Optimized performance

%tsum = polyval(wrev(c),Eps);
%tsum = polyval(wrev(c(1:N+1)),Eps);
tsum = polyval(c(N+1:-1:1),Eps);
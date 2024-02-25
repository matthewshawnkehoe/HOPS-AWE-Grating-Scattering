function [tsum] = taylorsum(c,Eps,N)
% taylorsum - Sums a truncated Taylor series.
%
% Inputs:
%
% c - Taylor series coefficients: [c_0,...,c_N]
% Eps - Value at which to sum
% N - Degree of truncated Taylor series (unused!)
%
% Outputs:
%
% tsum - Taylor sum evaluated at Eps
%
% DPN 2/7/12

%tsum = polyval(wrev(c),Eps);
%tsum = polyval(wrev(c(1:N+1)),Eps);
tsum = polyval(c(N+1:-1:1),Eps);
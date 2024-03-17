function [tsum] = taylorsum2(c,Eps,delta,N,M)
% taylorsum2.m: Performs a double summation of a truncated Taylor series by 
% using the coefficients.
%
% Inputs:
%   c: Taylor series coefficients for [c_{0,0},...,c_{0,N},...,c_{M,0},...,c_{M,N}]
%   Eps,delta: Value at which to sum
%   N,M: Degree of truncated Taylor series
%
% Outputs:
%   tsum: Taylor sum evaluated at Eps, delta
%
% MSK 7/11/2021

vander = (Eps.^(0:N).').*(delta.^(0:M));
tsum = sum(sum(c.*vander));

% Alternative Implementation (slightly faster for large M,N)
% c1 = zeros(N+1,1);
% c2 = zeros(1,M+1);
% c1(1,1)=1;
% c2(1,1)=1;
% c1(2:end) = cumprod(repmat(Eps,N,1));
% c2(2:end) = cumprod(repmat(delta,1,M));
% tsum = sum(c.*c1.*c2, 'all');

return;
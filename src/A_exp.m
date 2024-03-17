function [A_nm] = A_exp(alpha_m,f_x,N,M)
% A_exp.m: Computes A_nm as a reference for the field solvers.
%
%  Inputs:
%   alpha_m: a numerical constant 
%   f_x: the derivative of a test function in the x component
%   N: the maximum number of Taylor orders for the interfacial perturbation
%   M: the maximum number of Taylor orders for the frequency perturbation
%
%  Output:
%   A_nm: a numerical function passed into the field solvers.

Nx = length(f_x);
A_nm = zeros(Nx,N+1,M+1);

% n = 0, m = 0
A_nm(:,0+1,0+1) = ones(Nx,1);

% n = 0, m > 0: Do nothing as A_nm(:,0+1,m+1) = 0

% n > 0, m >= 0

for n=1:N
  for m=0:M
    sum = zeros(Nx,1);
    for r=0:m
      sum = sum + A_nm(:,n-1+1,r+1)*(1i*alpha_m(m-r+1));
    end
    A_nm(:,n+1,m+1) = f_x.*sum/n;
  end
end

return;
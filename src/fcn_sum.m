function [f] = fcn_sum(SumType,f_n_m,Eps,delta,...
    Nx,N,M)
% fcn_sum.m: Calculates the Taylor or Padé function sum.
%
%  Inputs:
%   SumType: boolean to control Taylor or Padé summation
%   f_n_m: a tensor representing the function values
%   Eps: the physical error in the surface deformation
%   delta: the numerical error in the discretization of the frequency perturbation
%   Nx: the number of discretization points
%   N: the maximum number of Taylor orders for the interfacial perturbation
%   M: the maximum number of Taylor orders for the frequency perturbation
%
%  Outputs:
%   f: the calculated function sum

% MSK 7/26/21: Changed f_n_m(j,r+1,s+1) to f_n_m(j,s+1,r+1)
f = zeros(Nx,1);
coeff = zeros(N+1,M+1);

for j=1:Nx
  for r=0:N
    for s=0:M
      coeff(r+1,s+1) = f_n_m(j,s+1,r+1);
    end
  end
  if(SumType==1)
    % Change to taylorsum_2_coeff if M,N >= 6
    f(j) = taylorsum_2_coeff(coeff,Eps,delta,N,M);
  elseif(SumType==2)
    f(j) = padesum2(coeff,Eps,delta,N,M);
  elseif(SumType==3)
    f(j) = padesum2_safe(coeff,Eps,delta,N,M);
  end
end

return;
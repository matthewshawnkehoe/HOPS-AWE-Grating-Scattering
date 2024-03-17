function [f] = fcn_sum_fast(SumType,f_n_m,Eps,delta,Nx,N,M)
% fcn_sum_fast.m: Efficiently calculates the Taylor or Padé function sum.
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

f = zeros(Nx,1);

for j=1:Nx
  if(SumType==1)
    % Change to taylorsum_2_coeff if M,N >= 6
    f(j) = taylorsum_2_coeff(f_n_m(:,:,j),Eps,delta,N,M);
  elseif(SumType==2)
    f(j) = padesum2(f_n_m(:,:,j),Eps,delta,N,M);
  else 
    f(j) = padesum2_safe(f_n_m(:,:,j),Eps,delta,N,M);
  end
end

return;
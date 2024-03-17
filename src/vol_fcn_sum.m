function [u] = vol_fcn_sum(SumType,u_n_m,Eps,delta,Nx,Nz,N,M)
% vol_fcn_sum.m: Calculates the Taylor or Padé volume function sum.
%
%  Inputs:
%   SumType: boolean to control Taylor or Padé summation
%   u_n_m: a tensor representing the approximate solution in the upper
%   field
%   Eps: the physical error in the surface deformation
%   delta: the numerical error in the discretization of the frequency perturbation
%   Nx: the number of discretization points
%   Nz: the number of collocation points
%   N: the maximum number of Taylor orders for the interfacial perturbation
%   M: the maximum number of Taylor orders for the frequency perturbation
%
%  Output:
%   u: the calculated volume function sum

u = zeros(Nx,Nz+1);
coeff = zeros(N+1,M+1);

% MSK 7/26/21: Changed u_n_m(j,ell+1,r+1,s+1) to u_n_m(j,ell+1,s+1,r+1)

for j=1:Nx
  for ell=0:Nz
    for r=0:N
      for s=0:M
        coeff(r+1,s+1) = u_n_m(j,ell+1,s+1,r+1);
      end
    end
    if(SumType==1)
      u(j,ell+1) = taylorsum2(coeff,Eps,delta,N,M);
    elseif(SumType==2)
      u(j,ell+1) = padesum2(coeff,Eps,delta,N,M);
    elseif(SumType==3)
      u(j,ell+1) = padesum2_safe(coeff,Eps,delta,N,M);
    end
  end
end

return;
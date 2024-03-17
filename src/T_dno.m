function [T] = T_dno(alpha,alphap,gamma,gammap,k2,Nx,M)
% T_dno.m: Calculates the frequency expansions of gamma_p per the perturbation in delta. 
% See Diffraction Problems: A HOPS/AWE Method Section 4.1, (12) - (16).
%
%  Inputs:
%   alpha: a numerical parameter
%   alphap: a numerical parameter calculated at every wave number p
%   gammap: a numerical parameter calculated at every wave number p
%   k2: alphap(0+1)^2 + gammap(0+1)^2
%   Nx: the number of discretization points 
%   M: the maximum number of Taylor orders for the frequency perturbation
%
%  Output:
%   T: vector resulting from the Taylor expansion of gamma_p

gammap_m = zeros(Nx,M+1);
gammap_m(:,0+1) = gammap;
gammap_m(:,1+1) = 2*(k2 - alpha*alphap)./(2*gammap);
gammap_m(:,2+1) = (gamma^2 - gammap_m(:,1+1).^2)./(2*gammap);

for m=3:M
  for r=1:m-1
    gammap_m(:,m+1) = gammap_m(:,m+1) ...
        - (gammap_m(:,m-r+1).*gammap_m(:,r+1))./(2*gammap);
  end
end

T = 1i*gammap_m;

return;
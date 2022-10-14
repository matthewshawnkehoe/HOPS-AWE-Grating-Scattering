function [T] = T_dno(alpha,alphap,gamma,gammap,k2,Nx,M)
% Tw_dno - Calculates the frequency expansions of gamma_p per the
% perbutation in delta. See Diffraction Problems: A HOPS/AWE Method
% section 4.1, (12) - (16).
%
% Inputs:
%
% alpha - parameter
% alphap - alpha + tilde(p) = alpha + (2pi / d) * p
% gamma - parameter
% gammap - sqrt(k^2 - alpha_p^2) provided p in propogating nodes
% k2 - alphap(0+1)^2 + gammap(0+1)^2;
% Nx - number of grid points evaluated on x-axis
% M - number of iterations for delta pertubation
%
% Outputs:
%
% Tw - Result of taylor expansion of gamma_p
%
% MSK 4/22/20

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
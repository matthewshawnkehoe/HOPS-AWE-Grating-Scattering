function [zeta_n_m,psi_n_m] = setup_zeta_psi_n_m(xx,pp,...
    alpha_bar,gamma_u_bar,f,f_x,Nx,N,M)
% setup_zeta_psi_n_m.m: Setup variables for calculating the difference in 
% the approximate solution (zeta_n_m) and the difference of the derivatives
% of the approximate solution (psi_n_m) in the upper and lower layers for
% validating the two layer solver.
%
%  Inputs:
%   xx: numerical discretization based on the number of discretization points in Nx
%   pp: an integer where tilde_p = (2*pi/d)*p and d is the periodicity of the grating interface
%   alpha_bar: a numerical constant
%   gamma_bar: a numerical constant
%   f: a test function representing the grating surface
%   f_x: the derivative of a test function in the x component
%   Nx: the number of discretization points
%   N: the maximum number of Taylor orders for the interfacial perturbation
%   M: the maximum number of Taylor orders for the frequency perturbation
%
%  Outputs:
%   zeta_n_m: tensor representing the approximate solution of the difference
%             between xi_u and xi_w at the surface z=g(x)
%   psi_n_m: tensor representing the approximate solution of the difference
%             between nu_u and nu_w at the surface z=g(x)

% MSK 7/26/21: Changed the size of zeta_n_m,psi_n_m from (Nx,N+1,M+1) to (Nx,M+1,N+1)

gamma_u_m = zeros(1,M+1);
gamma_u_m(0+1) = gamma_u_bar;
gamma_u_m(1+1) = gamma_u_bar;

alpha_u_m = zeros(1,M+1);
alpha_u_m(0+1) = alpha_bar;
alpha_u_m(1+1) = alpha_bar;

E_n_m = E_exp_lf(gamma_u_m,f,N,M);

zeta_n_m = zeros(Nx,M+1,N+1);
psi_n_m = zeros(Nx,M+1,N+1);

for n=0:N
  for m=0:M
    zeta_n_m(:,m+1,n+1) = -E_n_m(:,m+1,n+1);
    for ell=0:m
      psi_n_m(:,m+1,n+1) = psi_n_m(:,m+1,n+1) ...
          + (1i*gamma_u_m(m-ell+1)).*E_n_m(:,ell+1,n+1);
    end
    if(n>0)
      psi_n_m(:,m+1,n+1) = psi_n_m(:,m+1,n+1)...
          + f_x.*(1i*alpha_u_m(m-ell+1)).*E_n_m(:,ell+1,n-1+1);
    end    
  end
end

return;
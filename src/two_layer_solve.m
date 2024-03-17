function [U_n_m,W_n_m,ubar_n_m,wbar_n_m] ...
    = two_layer_solve(tau2,zeta_n_m,psi_n_m,...
    gamma_u_bar_p,gamma_w_bar_p,N,Nx,f,f_x,pp,alpha_bar,...
    gamma_u_bar,gamma_w_bar,Dz,a,b,Nz,M,identy,alpha_bar_p)
% two_layer_solve.m: Solves the two-layer problem in the upper and lower fields.
%
%  Inputs:
%   tau2: a numerical constant representing TE or TM mode
%   zeta_n_m: a tensor representing the Dirichlet data
%   psi_n_m: a tensor representing the Neumannn data
%   gamma_u_bar_p: a numerical parameter calculated at every wave number p in the upper field
%   gamma_w_bar_p: a numerical parameter calculated at every wave number p in the lower field
%   N: the maximum number of Taylor orders for the interfacial perturbation
%   Nx: the number of discretization points 
%   f: a test function representing the grating surface
%   f_x: the derivative of the test function in the x component
%   pp: an integer where tilde_p = (2*pi/d)*p and d is the periodicity of the grating interface
%   alpha_bar: a numerical constant
%   gamma_u_bar: a numerical constant in the upper field 
%   gamma_w_bar: a numerical constant in the lower field 
%   Dz: the partial derivative with respect to the z component
%   a: the artificial boundary imposed at the top of the upper layer
%   b: the artificial boundary imposed at the bottom of the lower layer
%   Nz: the number of collocation points
%   M: the maximum number of Taylor orders for the frequency perturbation
%   identy: the identity matrix
%   alpha_bar_p: a numerical constant at all wave numbers p
%
%  Outputs:
%   U_n_m: Approximate solution at the surface z=g(x) in the upper field
%   W_n_m: Approximate solution at the surface z=g(x) in the lower field
%   ubar_n_m: Approximate solution from the upper field solver at the top
%             collocation point
%   wbar_n_m: Approximate solution from the lower field solver at the lower
%             collocation point

% MSK 7/26/21: Changed the size of U_n_m,W_n_m,ubar_n_m, and wbar_n_m
% from (Nx,N+1,M+1) to (Nx,M+1,N+1)

U_n_m = zeros(Nx,M+1,N+1);
W_n_m = zeros(Nx,M+1,N+1);
ubar_n_m = zeros(Nx,M+1,N+1);
wbar_n_m = zeros(Nx,M+1,N+1);
                       
% n=0, m=0

Q = zeta_n_m(:,0+1,0+1);
R = -psi_n_m(:,0+1,0+1);
[U_n_m(:,0+1,0+1),W_n_m(:,0+1,0+1)] ...
    = AInverse(Q,R,gamma_u_bar_p,gamma_w_bar_p,Nx,tau2);
  
for n=0:N
  for m=0:M
    Q = zeta_n_m(:,m+1,n+1);
    R = -psi_n_m(:,m+1,n+1);
    % MSK 11/01/2021 - Add terms to R from removing the phase
    if n > 1 && m > 0
      R = R - f_x*1i*alpha_bar.*U_n_m(:,m,n-1) + tau2*f_x*1i*alpha_bar.*W_n_m(:,m,n-1);
      if m > 1
        R = R - f_x*1i*alpha_bar.*U_n_m(:,m-1,n-1) + tau2*f_x*1i*alpha_bar.*W_n_m(:,m-1,n-1);
      end
    end
    for r=0:n
      for s=0:m
        if((r<n)||(s<m))
          xi_n_m = zeros(Nx,M+1,N+1);
          xi_n_m(:,0+1,0+1) = U_n_m(:,s+1,r+1);
          u_n_m = field_tfe_helmholtz_m_and_n(xi_n_m,f,pp,gamma_u_bar_p,...
              alpha_bar,gamma_u_bar,Dz,a,Nx,Nz,n-r,m-s,identy,alpha_bar_p);
          G_n_m = dno_tfe_helmholtz_m_and_n(u_n_m,...
            f,pp,Dz,a,Nx,Nz,n-r,m-s);
          R = R - G_n_m(:,m-s+1,n-r+1);
        
          xi_n_m = zeros(Nx,M+1,N+1);
          xi_n_m(:,0+1,0+1) = W_n_m(:,s+1,r+1);
          w_n_m = field_tfe_helmholtz_m_and_n_lf(xi_n_m,f,pp,gamma_w_bar_p,...
              alpha_bar,gamma_w_bar,Dz,b,Nx,Nz,n-r,m-s,identy,alpha_bar_p);
          J_n_m = dno_tfe_helmholtz_m_and_n_lf(w_n_m,...
            f,pp,Dz,b,Nx,Nz,n-r,m-s);
          R = R - tau2*J_n_m(:,m-s+1,n-r+1);
        end
      end
    end
    if((n>0)||(m>0))
      [U_n_m(:,m+1,n+1),W_n_m(:,m+1,n+1)] ...
        = AInverse(Q,R,gamma_u_bar_p,gamma_w_bar_p,Nx,tau2);
    end
  end
end

u_n_m = field_tfe_helmholtz_m_and_n(U_n_m,f,pp,gamma_u_bar_p,...
    alpha_bar,gamma_u_bar,Dz,a,Nx,Nz,N,M,identy,alpha_bar_p);
w_n_m = field_tfe_helmholtz_m_and_n_lf(W_n_m,f,pp,gamma_w_bar_p,...
    alpha_bar,gamma_w_bar,Dz,b,Nx,Nz,N,M,identy,alpha_bar_p);

ell_top = 0 + 1;
ell_bottom = Nz + 1;

for n=0:N
  for m=0:M
    ubar_n_m(:,m+1,n+1) = u_n_m(:,ell_top,m+1,n+1);
    wbar_n_m(:,m+1,n+1) = w_n_m(:,ell_bottom,m+1,n+1);
  end
end

return;
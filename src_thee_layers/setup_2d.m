function [xx,kk,alphap,betap,eep,eem] = setup_2d(Nx,L,alpha,beta)
% setup_2d.m: Calculates the proper angle values and discretization for the
% numerical constants.
%
%  Inputs:
%   Nx: the number of discretization points 
%   L: the periodicity of the grating surface
%   alpha: a numerical constant
%   beta: a numerical constant
%
%  Outputs:
%   xx: numerical discretization based on the number of discretization points in Nx
%   kk: numerical discretization based on the number of discretization points in Nx
%   alphap: numerical constant based on discretization scheme
%   betap: numerical constant based on discretization scheme
%   eep: numerical constant representing a deviation from in the incidence radiation
%   eem: numerical constant representing a deviation from in the incidence radiation

xx = (L/Nx)*[0:Nx-1]';
kk = (2.0*pi/L)*[0:Nx/2-1,-Nx/2:-1]';
  
alphap = alpha + kk;
betap = 0*alphap;

kappa = sqrt(alpha^2 + beta^2);

value = kappa^2 - alphap.^2;
rho = abs(value);
theta = angle(value);

% 'angle' returns -pi < theta < pi
% If theta<0 add 2*pi so that 0 < theta < 2*pi

for j=1:Nx
  if(theta(j)<0)
    theta(j) = theta(j) + 2*pi;
  end
end

% Always choose betap so that Im{betap}>0
% Thus exp(i*betap*y) bounded for y->infty
% Thus exp(-i*betap*y) bounded for y->-infty

for j=1:Nx
  betap(j) = sqrt(rho(j))*exp(1i*theta(j)/2.0);
end

eep = exp(1i*alpha*xx);
eem = exp(-1i*alpha*xx);
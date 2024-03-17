function [U,W] = AInverse(Q,R,gammap,gammapw,Nx,tau2)
% AInverse.m: Computes the inverse of the matrix A.
%
%  Inputs:
%   Q: a numerical approximation of zeta_n_m 
%   R: a numerical approximation of psi_n_m 
%   gammap: a numerical constant in the upper field
%   gammapw: a numerical constant in the lower field
%   Nx: the number of discretization points 
%   tau2: a numerical constant representing TE or TM mode
%
%  Outputs:
%   U: a tensor representing the solution in the upper field at the interface
%   W: a tensor representing the solution in the lower field at the interface
  
a = zeros(Nx,1);
b = zeros(Nx,1);

Q_hat = fft(Q);
R_hat = fft(R);

for j=1:Nx
  det_p = tau2*gammapw(j) + gammap(j);
  a(j) = ((tau2*gammapw(j))*Q_hat(j) + 1i*R_hat(j))/det_p;
  b(j) = ((-gammap(j))*Q_hat(j) + 1i*R_hat(j))/det_p;
end

U = ifft(a);
W = ifft(b);

return;
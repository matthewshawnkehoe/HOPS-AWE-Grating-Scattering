function [What] = solvebvp_colloc_fast_lf(b,alpha,beta,gamma,d_min,n_min,...
    r_min,d_max,n_max,r_max,Nx,identy,D,D2,D_start,D_end)
% solvebvp_colloc_fast_lf.m: Solves the two-point boundary value problem by 
% accelerating the classical computation of Ax=b in the lower field.
%
%  Inputs:
%   b: the artificial boundary imposed at the bottom of the lower layer
%   alpha: a numerical constant
%   beta: a numerical constant
%   gamma: a numerical constant
%   d_min: the minimum in the two-point BVP for the d component 
%   n_min: the minimum in the two-point BVP for the n component 
%   r_min: the minimum in the two-point BVP for the r component 
%   d_max: the maximum in the two-point BVP for the d component 
%   n_max: the maximum in the two-point BVP for the n component
%   r_max: the maximum in the two-point BVP for the r component
%   Nx: the number of discretization points
%   identy: the identity matrix
%   D: rescaled Chebyshev differentiation matrix in computational domain
%   D2: square of the matrix D
%   D_start: start of the matrix D
%   D_end: end of the matrix D
%
%  Output:
%   W_hat: Fourier transform of the approximate solution in the lower field

% MSK 7/30/21: Created to vectorize most of the j loop

A = alpha*D2 + beta*D + reshape(gamma,1,1,Nx).*identy;
A(end,:,:) = repmat(n_min*D_end,[1,1,Nx]);
b(end,:) = r_min;
      
A(1,:,:) = repmat(n_max*D_start,[1,1,Nx]);
A(end,end,:) = A(end,end,:) + reshape(d_min,1,1,Nx);
A(1,1,:) = A(1,1,:) + d_max;
b(1,:) = r_max;

wtilde = pagemldivide(A,reshape(b,[],1,Nx));
What = wtilde(:,:).';

return;
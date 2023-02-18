function [Uhat] = solvebvp_colloc_fast(Uhat,b,alpha,beta,gamma,d_min,n_min,...
    r_min,d_max,n_max,r_max,Nx,identy,D,D2,D_start,D_end)

% MSK 7/30/21: Created to vectorize most of the j loop

A = alpha*D2 + beta*D + reshape(gamma,1,1,Nx).*identy;
A(end,:,:) = repmat(n_min*D_end,[1,1,Nx]);
b(end,:) = r_min;
      
A(1,:,:) = repmat(n_max*D_start,[1,1,Nx]);
A(end,end,:) = A(end,end,:) + d_min;
A(1,1,:) = A(1,1,:) + reshape(d_max,1,1,Nx);
b(1,:) = r_max;

utilde = pagemldivide(A,reshape(b,[],1,Nx));
Uhat = utilde(:,:).';

return;

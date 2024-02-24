function [utilde] = solvebvp_colloc(ftilde,alpha,beta,gamma,D,...
    d_a,n_a,r_a,d_b,n_b,r_b)

Ny = length(ftilde)-1;

D2 = D*D;
A = alpha*D2 + beta*D + gamma*eye(Ny+1);
b = ftilde;

A(Ny+1,:) = n_a*D(Ny+1,:);
A(Ny+1,Ny+1) = A(Ny+1,Ny+1) + d_a;
b(Ny+1) = r_a;

A(1,:) = n_b*D(1,:);
A(1,1) = A(1,1) + d_b;
b(1) = r_b;

utilde = A\b;

return;
function [Qnm] = iio_tfe_helmholtz_lower(unm, eta, f, p, alphap, gammap, ...
    eep, eem, Dz, a, Nx, Nz, N, M)
% iio_tfe_helmholtz_lower: Computes the intermediate output in the lower region for the Helmholtz equation.
%
%  Inputs:
%   unm: a tensor representing the solution at the previous step in the lower region
%   eta: a numerical constant
%   f: a test function representing the boundary
%   p: an integer where tilde_p = (2*pi/d)*p and d is the periodicity of the grating interface
%   alphap: a numerical constant at all wave numbers p
%   gammap: a numerical constant in the lower field for all wave numbers p
%   eep: a numerical constant
%   eem: a numerical constant
%   Dz: the partial derivative with respect to the z component
%   a: the artificial boundary imposed at the bottom of the lower layer
%   Nx: the number of discretization points
%   Nz: the number of collocation points
%   N: the maximum number of Taylor orders for the interfacial perturbation
%   M: the maximum number of Taylor orders for the frequency perturbation
%
%  Outputs:
%   Qnm: a tensor representing the intermediate output in the lower region

% Initialize the output tensor
Qnm = zeros(Nx, M+1, N+1);

ell_top = 1;
f_x = ifft((1i*p) .* fft(f));

for n = 0:N
    for m = 0:M
        u_z = dz(unm(:,:,m+1,n+1), Dz, a);
        Qnm(:,m+1,n+1) = u_z(:,ell_top) + 1i*eta*unm(:,ell_top,m+1,n+1);
        
        if n >= 1
            u_x = dxp(unm(:,:,m+1,n-1+1), alphap, eep, eem, Nx, Nz);
            Qnm(:,m+1,n+1) = Qnm(:,m+1,n+1) + f_x .* u_x(:,ell_top);
            Qnm(:,m+1,n+1) = Qnm(:,m+1,n+1) + (1i*eta/a) * f .* unm(:,ell_top,m+1,n-1+1);
            Qnm(:,m+1,n+1) = Qnm(:,m+1,n+1) + (1.0/a) * (f .* Qnm(:,m+1,n-1+1));
        end
        
        if n >= 2
            u_x = dxp(unm(:,:,m+1,n-2+1), alphap, eep, eem, Nx, Nz);
            Qnm(:,m+1,n+1) = Qnm(:,m+1,n+1) + (1.0/a) * (f .* (f_x .* u_x(:,ell_top)));
            u_z = dz(unm(:,:,m+1,n-1), Dz, a);
            Qnm(:,m+1,n+1) = Qnm(:,m+1,n+1) + f_x .* (f_x .* u_z(:,ell_top));
        end
    end
end

return;

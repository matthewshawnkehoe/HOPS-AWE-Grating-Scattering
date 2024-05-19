function [Runm, Rellnm] = iio_tfe_helmholtz_middle(unm, hbar, eta, fu, fell, ...
    p, alphap, gammap, eep, eem, Dz, a, Nx, Nz, N, M)
% iio_tfe_helmholtz_middle: Computes intermediate outputs in the Helmholtz equation using Taylor series expansions.
%
%  Inputs:
%   unm: a tensor representing the solution at the previous step
%   hbar: a numerical constant representing half the grid spacing in the z direction
%   eta: a numerical constant
%   fu: a test function representing the upper boundary
%   fell: a test function representing the lower boundary
%   p: an integer where tilde_p = (2*pi/d)*p and d is the periodicity of the grating interface
%   alphap: a numerical constant at all wave numbers p
%   gammap: a numerical constant in the upper field for all wave numbers p
%   eep: a numerical constant
%   eem: a numerical constant
%   Dz: the partial derivative with respect to the z component
%   a: the artificial boundary imposed at the top of the upper layer
%   Nx: the number of discretization points
%   Nz: the number of collocation points
%   N: the maximum number of Taylor orders for the interfacial perturbation
%   M: the maximum number of Taylor orders for the frequency perturbation
%
%  Outputs:
%   Runm: a tensor representing the intermediate solution at the upper boundary
%   Rellnm: a tensor representing the intermediate solution at the lower boundary

% Initialize the output tensors
Runm = zeros(Nx, M+1, N+1);
Rellnm = zeros(Nx, M+1, N+1);

ell_top = 1;
ell_bottom = Nz + 1;
fu_x = ifft((1i*p) .* fft(fu));
fell_x = ifft((1i*p) .* fft(fell));

for n = 0:N
    for m=0:M
        u_z = dz(unm(:,:,m+1,n+1), Dz, 2*hbar);
        Runm(:,m+1,n+1) = u_z(:,ell_top) + 1i*eta*unm(:,ell_top,m+1,n+1);
        Rellnm(:,m+1,n+1) = -u_z(:,ell_bottom) + 1i*eta*unm(:,ell_bottom,m+1,n+1);
        
        if n >= 1
            u_x = dxp(unm(:,:,m+1,n-1+1), alphap, eep, eem, Nx, Nz);
            
            Runm(:,m+1,n+1) = Runm(:,m+1,n+1) - fu_x .* u_x(:,ell_top);
            Runm(:,m+1,n+1) = Runm(:,m+1,n+1) + (1i*eta/(2*hbar)) * fu .* unm(:,ell_top,m+1,n-1+1);
            Runm(:,m+1,n+1) = Runm(:,m+1,n+1) - (1i*eta/(2*hbar)) * fell .* unm(:,ell_top,m+1,n-1+1);
            Runm(:,m+1,n+1) = Runm(:,m+1,n+1) - (1.0/(2*hbar)) * fu .* Runm(:,m+1,n+1);
            Runm(:,m+1,n+1) = Runm(:,m+1,n+1) + (1.0/(2*hbar)) * fell .* Runm(:,m+1,n+1);
            
            Rellnm(:,m+1,n+1) = Rellnm(:,m+1,n+1) + fell_x .* u_x(:,ell_bottom);
            Rellnm(:,m+1,n+1) = Rellnm(:,m+1,n+1) + (1i*eta/(2*hbar)) * fu .* unm(:,ell_bottom,m+1,n-1+1);
            Rellnm(:,m+1,n+1) = Rellnm(:,m+1,n+1) - (1i*eta/(2*hbar)) * fell .* unm(:,ell_bottom,m+1,n-1+1);
            Rellnm(:,m+1,n+1) = Rellnm(:,m+1,n+1) - (1.0/(2*hbar)) * fu .* Rellnm(:,m+1,n+1);
            Rellnm(:,m+1,n+1) = Rellnm(:,m+1,n+1) + (1.0/(2*hbar)) * fell .* Rellnm(:,m+1,n+1);
        end
        
        if n >= 2
            u_x = dxp(unm(:,:,m+1,n-2+1), alphap, eep, eem, Nx, Nz);
            u_z = dz(unm(:,:,m+1,n-2+1), Dz, 2*hbar);
            
            Runm(:,m+1,n+1) = Runm(:,m+1,n+1) - (1.0/(2*hbar)) * fu .* (fu_x .* u_x(:,ell_top));
            Runm(:,m+1,n+1) = Runm(:,m+1,n+1) + (1.0/(2*hbar)) * fell .* (fu_x .* u_x(:,ell_top));
            Runm(:,m+1,n+1) = Runm(:,m+1,n+1) + fu_x .* (fu_x .* u_z(:,ell_top));
            
            Rellnm(:,m+1,n+1) = Rellnm(:,m+1,n+1) + (1.0/(2*hbar)) * fu .* (fell_x .* u_x(:,ell_bottom));
            Rellnm(:,m+1,n+1) = Rellnm(:,m+1,n+1) - (1.0/(2*hbar)) * fell .* (fell_x .* u_x(:,ell_bottom));
            Rellnm(:,m+1,n+1) = Rellnm(:,m+1,n+1) - fell_x .* (fell_x .* u_z(:,ell_bottom));
        end
    end
end

return;

function [f,f_x] = fourier_repr_lipschitz(P,x)
% fourier_repr_lipschitz.m: Plots the Fourier representation of the 
% Lipschitz boundary profile as in Section 6.4 of "A Stable HOPS/AWE Method
% for Grating Scattering".
%
%  Inputs:
%   P: An integer representing the maximum summation order of a truncated
%      Fourier series (to minimize the effect of aliasing errors)
%   x: numerical discretization based on the number of discretization points in Nx
%
%  Outputs:
%   f: The Fourier representation of a test function at the grating surface
%   f_x: The Fourier representation of the partial derivative in the x
%        component for the test function at the grating surface

f = 0;
for k = 1:ceil(P/2)
    f = f + 8/(pi^2*(2*k-1)^2) * cos((2*k-1)*x);
end

f_x = 0;
for k = 1:ceil(P/2)
    f_x = f_x + (-8*(2*k-1))/(pi^2*(2*k-1)^2) * sin((2*k-1)*x);
end

figure(3);
hold on;
subplot(1,2,1)
plot(x,f);
title('Lipschitz profile for $f$','interpreter','latex','FontSize',16);
subplot(1,2,2)
plot(x,f_x);
title('Lipschitz profile for $f_x$','interpreter','latex','FontSize',16);

return
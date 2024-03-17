function [f,f_x] = fourier_repr_rough(P,x)
% fourier_repr_rough.m: Plots the Fourier representation of the moderately
% smooth (C^4) boundary profile as in Section 6.4 of "A Stable HOPS/AWE Method
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
for k = 1:P
    f = f + 96*(2*k^2*pi^2 - 21)/(125*k^8) * cos(k*x);
end

f_x = 0;
for k = 1:P
    f_x = f_x + (-96*k*(2*k^2*pi^2 - 21))/(125*k^8) * sin(k*x);
end

figure(4);
hold on;
subplot(1,2,1)
plot(x,f);
title('Rough profile for $f$','interpreter','latex','FontSize',16);
subplot(1,2,2)
plot(x,f_x);
title('Rough profile for $f_x$','interpreter','latex','FontSize',16);

return
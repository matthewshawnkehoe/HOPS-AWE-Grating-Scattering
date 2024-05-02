function [u_x] = dx(u,p)
% dx.m: Computes the partial derivative in the x component through the FFT and IFFT.
%
%  Inputs:
%   u: a tensor representing the approximate solution in a given layer
%   p: an integer where tilde_p = (2*pi/d)*p and d is the periodicity of the grating interface
%
%  Output:
%   u_x: a tensor representing the derivative of the approximate solution in the x component

u_x = ifft((1i*p).*fft(u));

return;
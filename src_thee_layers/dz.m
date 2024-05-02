function [u_z] = dz(u,Dz,b)
% dz.m: Computes the partial derivative in the z component.
%
%  Inputs:
%   u: a tensor representing the approximate solution in a given layer
%   p: an integer where tilde_p = (2*pi/d)*p and d is the periodicity of the grating interface
%   b: the artificial boundary imposed at the bottom of the lower layer
%
%  Output:
%   u_z: a tensor representing the derivative of the approximate solution in the z component

u_z = ((2.0/b)*Dz*u.').';

return;
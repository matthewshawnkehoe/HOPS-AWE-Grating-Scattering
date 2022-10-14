function [u_z] = dz(u,Dz,b)

u_z = ((2.0/b)*Dz*u.').';

return;
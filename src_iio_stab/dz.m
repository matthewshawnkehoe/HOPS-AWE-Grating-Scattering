function [u_z] = dz(u,Dz,b,Nx,Nz)

u_z = zeros(Nx,Nz+1);
g = zeros(Nz+1,1);

for j=1:Nx
  for ell=0:Nz
    g(ell+1) = u(j,ell+1);
  end
  u_z(j,:) = (2.0/b)*Dz*g;
end

return;
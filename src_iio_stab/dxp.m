function [u_x] = dxp(u,alphap,eep,eem,Nx,Ny)

u_x = zeros(Nx,Ny+1);

for ell=0:Ny
  f = u(:,ell+1);
  u_x(:,ell+1) = eep.*ifft( (1i*alphap).*fft(eem.*f) );
end

return;
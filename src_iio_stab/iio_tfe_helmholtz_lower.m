function [Sn] = iio_tfe_helmholtz_lower(wn,eta,f,p,alphap,gammap,...
    eep,eem,Dz,b,Nx,Nz,N)

Sn = zeros(Nx,N+1);

ell_top = 0 + 1;
f_x = ifft( (1i*p).*fft(f) );

for n=0:N
  w_z = dz(wn(:,:,n+1),Dz,b,Nx,Nz);
  Sn(:,n+1) = w_z(:,ell_top) + 1i*eta*wn(:,ell_top,n+1);
  if(n>=1)
    w_x = dxp(wn(:,:,n-1+1),alphap,eep,eem,Nx,Nz);
    Sn(:,n+1) = Sn(:,n+1) - f_x.*w_x(:,ell_top);
    
    Sn(:,n+1) = Sn(:,n+1) + (1i*eta/b)*f.*wn(:,ell_top,n-1+1);
    
    Sn(:,n+1) = Sn(:,n+1) - (1.0/b)*(f.*Sn(:,n-1+1));
  end
  if(n>=2)
    w_x = dxp(wn(:,:,n-2+1),alphap,eep,eem,Nx,Nz);
    Sn(:,n+1) = Sn(:,n+1) - (1.0/b)*(f.*(f_x.*w_x(:,ell_top)));

    w_z = dz(wn(:,:,n-2+1),Dz,b,Nx,Nz);
    Sn(:,n+1) = Sn(:,n+1) + f_x.*(f_x.*w_z(:,ell_top));
  end
end

return;
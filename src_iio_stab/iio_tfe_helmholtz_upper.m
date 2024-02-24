function [Qn] = iio_tfe_helmholtz_upper(un,eta,f,p,alphap,gammap,...
    eep,eem,Dz,a,Nx,Nz,N)

Qn = zeros(Nx,N+1);

ell_bottom = Nz + 1;
f_x = ifft( (1i*p).*fft(f) );

for n=0:N
  u_z = dz(un(:,:,n+1),Dz,a,Nx,Nz);
  Qn(:,n+1) = -u_z(:,ell_bottom) + 1i*eta*un(:,ell_bottom,n+1);
  if(n>=1)
    u_x = dxp(un(:,:,n-1+1),alphap,eep,eem,Nx,Nz);
    Qn(:,n+1) = Qn(:,n+1) + f_x.*u_x(:,ell_bottom);
    
    Qn(:,n+1) = Qn(:,n+1) - (1i*eta/a)*f.*un(:,ell_bottom,n-1+1);
    
    Qn(:,n+1) = Qn(:,n+1) + (1.0/a)*(f.*Qn(:,n-1+1));
  end
  if(n>=2)
    u_x = dxp(un(:,:,n-2+1),alphap,eep,eem,Nx,Nz);
    Qn(:,n+1) = Qn(:,n+1) - (1.0/a)*(f.*(f_x.*u_x(:,ell_bottom)));

    u_z = dz(un(:,:,n-2+1),Dz,a,Nx,Nz);
    Qn(:,n+1) = Qn(:,n+1) - f_x.*(f_x.*u_z(:,ell_bottom));
  end
end

return;

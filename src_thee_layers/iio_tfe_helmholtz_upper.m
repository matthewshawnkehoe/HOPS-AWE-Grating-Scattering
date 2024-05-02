function [Qnm] = iio_tfe_helmholtz_upper(unm,eta,f,p,alphap,gammap,...
    eep,eem,Dz,a,Nx,Nz,N,M)

Qnm = zeros(Nx,M+1,N+1);

ell_bottom = Nz + 1;
f_x = ifft( (1i*p).*fft(f) );

for n=0:N
  for m=0:M
    u_z = dz(unm(:,:,m+1,n+1),Dz,a);
    Qnm(:,m+1,n+1) = -u_z(:,ell_bottom) + 1i*eta*unm(:,ell_bottom,m+1,n+1);
    if(n>=1)
      u_x = dxp(unm(:,:,m+1,n-1+1),alphap,eep,eem,Nx,Nz);
      Qnm(:,m+1,n+1) = Qnm(:,m+1,n+1) + f_x.*u_x(:,ell_bottom);
    
      Qnm(:,m+1,n+1) = Qnm(:,m+1,n+1) - (1i*eta/a)*f.*unm(:,ell_bottom,m+1,n-1+1);
    
      Qnm(:,m+1,n+1) = Qnm(:,m+1,n+1) + (1.0/a)*(f.*Qnm(:,m+1,n-1+1));
    end
    if(n>=2)
      u_x = dxp(unm(:,:,m+1,n-2+1),alphap,eep,eem,Nx,Nz);
      Qnm(:,m+1,n+1) = Qnm(:,m+1,n+1) - (1.0/a)*(f.*(f_x.*u_x(:,ell_bottom)));

      u_z = dz(unm(:,:,m+1,n-2+1),Dz,a);
      Qnm(:,m+1,n+1) = Qnm(:,m+1,n+1) - f_x.*(f_x.*u_z(:,ell_bottom));
    end
  end
end

return;

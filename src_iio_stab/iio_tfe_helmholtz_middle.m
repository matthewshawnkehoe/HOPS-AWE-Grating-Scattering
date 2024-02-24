function [Run,Relln] = iio_tfe_helmholtz_middle(un,hbar,eta,fu,fell,...
    p,alphap,gammap,eep,eem,Dz,a,Nx,Nz,N)

Run = zeros(Nx,N+1);
Relln = zeros(Nx,N+1);

ell_top = 0 + 1;
ell_bottom = Nz + 1;
fu_x = ifft( (1i*p).*fft(fu) );
fell_x = ifft( (1i*p).*fft(fell) );

for n=0:N
  u_z = dz(un(:,:,n+1),Dz,2*hbar,Nx,Nz);
  Run(:,n+1) = u_z(:,ell_top) + 1i*eta*un(:,ell_top,n+1);
  Relln(:,n+1) = -u_z(:,ell_bottom) + 1i*eta*un(:,ell_bottom,n+1);
  if(n>=1)
    u_x = dxp(un(:,:,n-1+1),alphap,eep,eem,Nx,Nz);
    
    Run(:,n+1) = Run(:,n+1) - fu_x.*u_x(:,ell_top);
    Run(:,n+1) = Run(:,n+1) + (1i*eta/(2*hbar))*fu.*un(:,ell_top,n-1+1);
    Run(:,n+1) = Run(:,n+1) - (1i*eta/(2*hbar))*fell.*un(:,ell_top,n-1+1);
    Run(:,n+1) = Run(:,n+1) - (1.0/(2*hbar))*fu.*Run(:,n-1+1);
    Run(:,n+1) = Run(:,n+1) + (1.0/(2*hbar))*fell.*Run(:,n-1+1);
    
    Relln(:,n+1) = Relln(:,n+1) + fell_x.*u_x(:,ell_bottom);
    Relln(:,n+1) = Relln(:,n+1) + (1i*eta/(2*hbar))*fu.*un(:,ell_bottom,n-1+1);
    Relln(:,n+1) = Relln(:,n+1) - (1i*eta/(2*hbar))*fell.*un(:,ell_bottom,n-1+1);
    Relln(:,n+1) = Relln(:,n+1) - (1.0/(2*hbar))*fu.*Relln(:,n-1+1);
    Relln(:,n+1) = Relln(:,n+1) + (1.0/(2*hbar))*fell.*Relln(:,n-1+1);
  end
  if(n>=2)
    u_x = dxp(un(:,:,n-2+1),alphap,eep,eem,Nx,Nz);
    u_z = dz(un(:,:,n-2+1),Dz,2*hbar,Nx,Nz);
    
    Run(:,n+1) = Run(:,n+1) - (1.0/(2*hbar))*fu.*(fu_x.*u_x(:,ell_top));
    Run(:,n+1) = Run(:,n+1) + (1.0/(2*hbar))*fell.*(fu_x.*u_x(:,ell_top));
    Run(:,n+1) = Run(:,n+1) + fu_x.*(fu_x.*u_z(:,ell_top));
    
    Relln(:,n+1) = Relln(:,n+1) + (1.0/(2*hbar))*fu.*(fell_x.*u_x(:,ell_bottom));
    Relln(:,n+1) = Relln(:,n+1) - (1.0/(2*hbar))*fell.*(fell_x.*u_x(:,ell_bottom));
    Relln(:,n+1) = Relln(:,n+1) - fell_x.*(fell_x.*u_z(:,ell_bottom));
  end
end

return;
function [Gnm] = dno_tfe_helmholtz_m_and_n_lf(wnm,f,p,Dz,...
    b,Nx,Nz,N,M)

Gnm = zeros(Nx,M+1,N+1);

ell_top = 0 + 1;
f_x = ifft( (1i*p).*fft(f) );

for n=0:N
  for m=0:M
    w_z = dz(wnm(:,:,m+1,n+1),Dz,b);
    Gnm(:,m+1,n+1) = w_z(:,ell_top);
    if(n>=1)
      w_x = dx(wnm(:,:,m+1,n-1+1),p);
      Gnm(:,m+1,n+1) = Gnm(:,m+1,n+1) - f_x.*w_x(:,ell_top);
    
      Gnm(:,m+1,n+1) = Gnm(:,m+1,n+1) - (1.0/b)*(f.*Gnm(:,m+1,n-1+1));
    end
    if(n>=2)
      w_x = dx(wnm(:,:,m+1,n-2+1),p);
      Gnm(:,m+1,n+1) = Gnm(:,m+1,n+1) - (1.0/b)*(f.*(f_x.*w_x(:,ell_top)));

      w_z = dz(wnm(:,:,m+1,n-2+1),Dz,b);
      Gnm(:,m+1,n+1) = Gnm(:,m+1,n+1) + f_x.*(f_x.*w_z(:,ell_top));
    end
  end
end

return;
function [wn] = field_tfe_helmholtz_lower(Wn,eta,f,p,alphap,gammap,...
    eep,eem,Dz,b,Nx,Nz,N)

wn = zeros(Nx,Nz+1,N+1);
Wnhat = zeros(Nx,N+1);

k2 = alphap(0+1)^2 + gammap(0+1)^2;

ell_top = 0 + 1;
ell_bottom = Nz + 1;
for n=0:N
  Wnhat(:,n+1) = fft(eem.*Wn(:,n+1));
end
f_x = real(ifft( (1i*p).*fft(f) ));

ll = [0:Nz]';
z_min = -b; z_max = 0.0;
z = ((z_max-z_min)/2.0)*(cos(pi*ll/Nz) - 1.0) + z_max;

f_full = zeros(Nx,Nz+1);
f_x_full = zeros(Nx,Nz+1);
for ell=0:Nz
  f_full(:,ell+1) = f;
  f_x_full(:,ell+1) = f_x;
end
temphat = zeros(Nx,Nz+1);

z_plus_b_full = zeros(Nx,Nz+1);
for j=1:Nx
  z_plus_b_full(j,:) = z.' + b;
end

% Order zero

A0 = -1i*gammap - 1i*eta;
for ell=0:Nz
  wn(:,ell+1,0+1) = eep.*ifft( exp(-1i*gammap*z(ell+1)).*Wnhat(:,0+1)./A0 );
end

% Order n>0

for n=1:N
    
  % Form Fn, Jn
  
  Fn = zeros(Nx,Nz+1);
  Hn = Wn(:,n+1);
  Jn = zeros(Nx,1);
  
  A1_xx = (2.0/b)*f_full;
  A1_xz = -(1.0/b)*(z_plus_b_full).*f_x_full;
  A1_zx = A1_xz;
  %A1_zz = 0;
  
  A2_xx = (1.0/b^2)*f_full.^2;
  A2_xz = -(1.0/b^2)*(z_plus_b_full).*(f_full.*f_x_full);
  A2_zx = A2_xz;
  A2_zz = (1.0/b^2)*((z_plus_b_full).^2).*(f_x_full.^2);
  
  B1_x = (1.0/b)*f_x_full;
  %B1_z = 0;
  
  B2_x = (1.0/b^2)*f_full.*f_x_full;
  B2_z = -(1.0/b^2).*(z_plus_b_full).*(f_x_full.^2);
  
  C1 = k2*(2.0/b)*f_full;
  C2 = k2*(1.0/b^2)*f_full.^2;
  
  if(n>=1)
    w_x = dxp(wn(:,:,n-1+1),alphap,eep,eem,Nx,Nz);
    temp = A1_xx.*w_x;
    Fn = Fn - dxp(temp,alphap,eep,eem,Nx,Nz);
    temp = A1_zx.*w_x;
    Fn = Fn - dz(temp,Dz,b,Nx,Nz);
    temp = B1_x.*w_x;
    Fn = Fn + temp;
    
    w_z = dz(wn(:,:,n-1+1),Dz,b,Nx,Nz);
    temp = A1_xz.*w_z;
    Fn = Fn - dxp(temp,alphap,eep,eem,Nx,Nz);
    %A1_zz = 0
    %B1_z = 0
    
    temp = C1.*wn(:,:,n-1+1);
    Fn = Fn - temp;
    
    Su = eep.*ifft( (-1i*gammap).*fft(eem.*wn(:,ell_bottom,n-1+1)) );
    Jn = Jn + (1.0/b)*f.*Su;
    Hn = Hn + (1.0/b)*f.*Wn(:,n-1+1)...
        + (1i*eta/b)*f.*wn(:,ell_top,n-1+1)...
        + f_x.*w_x(:,ell_top);
  end
  
  if(n>=2)
    w_x = dxp(wn(:,:,n-2+1),alphap,eep,eem,Nx,Nz);
    temp = A2_xx.*w_x;
    Fn = Fn - dxp(temp,alphap,eep,eem,Nx,Nz);
    temp = A2_zx.*w_x;
    Fn = Fn - dz(temp,Dz,b,Nx,Nz);
    temp = B2_x.*w_x;
    Fn = Fn + temp;
    
    w_z = dz(wn(:,:,n-2+1),Dz,b,Nx,Nz);
    temp = A2_xz.*w_z;
    Fn = Fn - dxp(temp,alphap,eep,eem,Nx,Nz);
    temp = A2_zz.*w_z;
    Fn = Fn - dz(temp,Dz,b,Nx,Nz);
    temp = B2_z.*w_z;
    Fn = Fn + temp;
    
    temp = C2.*wn(:,:,n-2+1);
    Fn = Fn - temp;
    
    Hn = Hn + (1.0/b)*f.*f_x.*w_x(:,ell_top)...
        - f_x.*f_x.*w_z(:,ell_top);
  end
  
  % Solve elliptic equation
  
  Fnhat = zeros(Nx,Nz+1);
  for ell=0:Nz
    Fnhat(:,ell+1) = fft(eem.*Fn(:,ell+1));
  end
  Hnhat = fft(eem.*Hn);
  Jnhat = fft(eem.*Jn);
  
  for j=1:Nx
    Fnhat_p = Fnhat(j,:).';
    alphaalpha = 1.0;
    betabeta = 0.0;
    gammagamma = k2 - (alphap(j))^2;
    d_a = -(-1i*gammap(j));
    n_a = 1.0;
    r_a = Jnhat(j);
    d_b = -1i*eta;
    n_b = 1.0;
    r_b = Hnhat(j);
    what_p = solvebvp_colloc(Fnhat_p,alphaalpha,betabeta,gammagamma,...
        (2.0/(z_max-z_min))*Dz,d_a,n_a,r_a,d_b,n_b,r_b);
    
    temphat(j,:) = what_p.';
  end
  
  for ell=0:Nz
    wn(:,ell+1,n+1) = eep.*ifft(temphat(:,ell+1));
  end

end

return;
function [un] = field_tfe_helmholtz_upper(Un,eta,f,p,alphap,gammap,...
    eep,eem,Dz,a,Nx,Nz,N)

un = zeros(Nx,Nz+1,N+1);
Unhat = zeros(Nx,N+1);

k2 = alphap(0+1)^2 + gammap(0+1)^2;

ell_top = 0 + 1;
ell_bottom = Nz + 1;
for n=0:N
  Unhat(:,n+1) = fft(eem.*Un(:,n+1));
end
f_x = real(ifft( (1i*p).*fft(f) ));

ll = [0:Nz]';
z_min = 0.0; z_max = a;
z = ((z_max-z_min)/2.0)*(cos(pi*ll/Nz) - 1.0) + z_max;

f_full = zeros(Nx,Nz+1);
f_x_full = zeros(Nx,Nz+1);
for ell=0:Nz
  f_full(:,ell+1) = f;
  f_x_full(:,ell+1) = f_x;
end
temphat = zeros(Nx,Nz+1);

a_minus_z_full = zeros(Nx,Nz+1);
for j=1:Nx
  a_minus_z_full(j,:) = a - z.';
end

% Order zero

A0 = -1i*gammap - 1i*eta;
for ell=0:Nz
  un(:,ell+1,0+1) = eep.*ifft( exp(1i*gammap*z(ell+1)).*Unhat(:,0+1)./A0 );
end

% Order n>0

for n=1:N
    
  % Form Fn, Jn
  
  Fn = zeros(Nx,Nz+1);
  Hn = Un(:,n+1);
  Jn = zeros(Nx,1);
  
  A1_xx = -(2.0/a)*f_full;
  A1_xz = -(1.0/a)*(a_minus_z_full).*f_x_full;
  A1_zx = A1_xz;
  %A1_zz = 0;
  
  A2_xx = (1.0/a^2)*f_full.^2;
  A2_xz = (1.0/a^2)*(a_minus_z_full).*(f_full.*f_x_full);
  A2_zx = A2_xz;
  A2_zz = (1.0/a^2)*((a_minus_z_full).^2).*(f_x_full.^2);
  
  B1_x = -(1.0/a)*f_x_full;
  %B1_z = 0;
  
  B2_x = (1.0/a^2)*f_full.*f_x_full;
  B2_z = (1.0/a^2).*(a_minus_z_full).*(f_x_full.^2);
  
  C1 = -k2*(2.0/a)*f_full;
  C2 = k2*(1.0/a^2)*f_full.^2;
  
  if(n>=1)
    u_x = dxp(un(:,:,n-1+1),alphap,eep,eem,Nx,Nz);
    temp = A1_xx.*u_x;
    Fn = Fn - dxp(temp,alphap,eep,eem,Nx,Nz);
    temp = A1_zx.*u_x;
    Fn = Fn - dz(temp,Dz,a,Nx,Nz);
    temp = B1_x.*u_x;
    Fn = Fn + temp;
    
    u_z = dz(un(:,:,n-1+1),Dz,a,Nx,Nz);
    temp = A1_xz.*u_z;
    Fn = Fn - dxp(temp,alphap,eep,eem,Nx,Nz);
    %A1_zz = 0
    %B1_z = 0
    
    temp = C1.*un(:,:,n-1+1);
    Fn = Fn - temp;
    
    Su = eep.*ifft( (1i*gammap).*fft(eem.*un(:,ell_top,n-1+1)) );
    Jn = Jn - (1.0/a)*f.*Su;
    Hn = Hn - (1.0/a)*f.*Un(:,n-1+1)...
        - (1i*eta/a)*f.*un(:,ell_bottom,n-1+1)...
        - f_x.*u_x(:,ell_bottom);
  end
  
  if(n>=2)
    u_x = dxp(un(:,:,n-2+1),alphap,eep,eem,Nx,Nz);
    temp = A2_xx.*u_x;
    Fn = Fn - dxp(temp,alphap,eep,eem,Nx,Nz);
    temp = A2_zx.*u_x;
    Fn = Fn - dz(temp,Dz,a,Nx,Nz);
    temp = B2_x.*u_x;
    Fn = Fn + temp;
    
    u_z = dz(un(:,:,n-2+1),Dz,a,Nx,Nz);
    temp = A2_xz.*u_z;
    Fn = Fn - dxp(temp,alphap,eep,eem,Nx,Nz);
    temp = A2_zz.*u_z;
    Fn = Fn - dz(temp,Dz,a,Nx,Nz);
    temp = B2_z.*u_z;
    Fn = Fn + temp;
    
    temp = C2.*un(:,:,n-2+1);
    Fn = Fn - temp;
    
    Hn = Hn + (1.0/a)*f.*f_x.*u_x(:,ell_bottom)...
        + f_x.*f_x.*u_z(:,ell_bottom);
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
    d_a = -1i*eta;
    n_a = -1.0;
    r_a = Hnhat(j);
    d_b = -1i*gammap(j);
    n_b = 1.0;
    r_b = Jnhat(j);
    uhat_p = solvebvp_colloc(Fnhat_p,alphaalpha,betabeta,gammagamma,...
        (2.0/(z_max-z_min))*Dz,d_a,n_a,r_a,d_b,n_b,r_b);
    temphat(j,:) = uhat_p.';
  end
  
  for ell=0:Nz
    un(:,ell+1,n+1) = eep.*ifft(temphat(:,ell+1));
  end

end

return;
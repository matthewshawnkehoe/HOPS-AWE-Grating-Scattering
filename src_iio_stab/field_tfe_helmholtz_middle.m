function [un] = field_tfe_helmholtz_middle(Vun,Velln,hbar,eta,fu,fell,...
    p,alphap,gammap,eep,eem,Dz,a,Nx,Nz,N)

un = zeros(Nx,Nz+1,N+1);
Vunhat = zeros(Nx,N+1);
Vellnhat = zeros(Nx,N+1);

k2 = alphap(0+1)^2 + gammap(0+1)^2;

ell_top = 0 + 1;
ell_bottom = Nz + 1;
for n=0:N
  Vunhat(:,n+1) = fft(eem.*Vun(:,n+1));
  Vellnhat(:,n+1) = fft(eem.*Velln(:,n+1));
end
fu_x = real(ifft( (1i*p).*fft(fu) ));
fell_x = real(ifft( (1i*p).*fft(fell) ));

ll = [0:Nz]';
z_min = -hbar; z_max = hbar;
z = ((z_max-z_min)/2.0)*(cos(pi*ll/Nz) - 1.0) + z_max;

fu_full = zeros(Nx,Nz+1);
fell_full = zeros(Nx,Nz+1);
fu_x_full = zeros(Nx,Nz+1);
fell_x_full = zeros(Nx,Nz+1);
for ell=0:Nz
  fu_full(:,ell+1) = fu;
  fell_full(:,ell+1) = fell;
  fu_x_full(:,ell+1) = fu_x;
  fell_x_full(:,ell+1) = fell_x;
end
temphat = zeros(Nx,Nz+1);

Z_L = zeros(Nx,Nz+1);
Z_U = zeros(Nx,Nz+1);
for ell=0:Nz
  for j=1:Nx
    Z_L(j,ell+1) = (z(ell+1) - z_min)/(2.0*hbar);
    Z_U(j,ell+1) = (z_max - z(ell+1))/(2.0*hbar);
  end
end

% Order zero

sh = sinh(1i*gammap*hbar)./(1i*gammap);
ch = cosh(1i*gammap*hbar);
aa = -(gammap.^2).*sh - 1i*eta*ch;
bb = ch - 1i*eta*sh;
Bp = (Vunhat(:,0+1) + Vellnhat(:,0+1))./(2.0*aa);
Cp = (Vunhat(:,0+1) - Vellnhat(:,0+1))./(2.0*bb);
for ell=0:Nz
  un(:,ell+1,0+1) = eep.*ifft( Bp.*cosh(1i*gammap*z(ell+1))...
      + Cp.*sinh(1i*gammap*z(ell+1))./(1i*gammap) );
end

% Order n>0

for n=1:N
    
  % Form Fn, Jn
  
  Fn = zeros(Nx,Nz+1);
  Hun = Vun(:,n+1);
  Helln = Velln(:,n+1);
  
  A10_xx = (-2.0/(2.0*hbar))*fell_full;
  A10_xz = -Z_U.*fell_x_full;
  A10_zx = A10_xz;
  %A10_zz = 0;
  
  A01_xx = (2.0/(2.0*hbar))*fu_full;
  A01_xz = -Z_L.*fu_x_full;
  A01_zx = A01_xz;
  %A01_zz = 0;
  
  A1_xx = A10_xx + A01_xx;
  A1_xz = A10_xz + A01_xz;
  A1_zx = A10_zx + A01_zx;
  %A1_zz = 0;
  
  A20_xx = (1.0/(2*hbar)^2)*fell_full.^2;
  A20_xz = (1.0/(2*hbar))*Z_U.*fell_full.*fell_x_full;
  A20_zx = A20_xz;
  A20_zz = (Z_U.^2).*(fell_x_full).^2;
  
  A11_xx = -(2.0/(2*hbar)^2)*fell_full.*fu_full;
  A11_xz = (1.0/(2*hbar))*(Z_L.*fell_full.*fu_x_full...
      -Z_U.*fu_full.*fell_x_full);
  A11_zx = A11_xz;
  A11_zz = 2*Z_U.*Z_L.*fell_x_full.*fu_x_full;
  
  A02_xx = (1.0/(2*hbar)^2)*fu_full.^2;
  A02_xz = -(1.0/(2*hbar))*Z_L.*fu_full.*fu_x_full;
  A02_zx = A02_xz;
  A02_zz = (Z_L.^2).*(fu_x_full).^2;
    
  A2_xx = A20_xx + A11_xx + A02_xx;
  A2_xz = A20_xz + A11_xz + A02_xz;
  A2_zx = A2_xz;
  A2_zz = A20_zz + A11_zz + A02_zz;
  
  B10_x = (1.0/(2.0*hbar))*fell_x_full;
  %B10_z = 0
  B01_x = -(1.0/(2.0*hbar))*fu_x_full;
  %B01_z = 0
  B20_x = -(1.0/(2.0*hbar)^2)*fell_full.*fell_x_full;
  B20_z = -(1.0/(2.0*hbar))*Z_U.*fell_x_full.^2;
  B11_x = (1.0/(2.0*hbar)^2)*(fu_full.*fell_x_full + fell_full.*fu_x_full);
  B11_z = (1.0/(2.0*hbar))*(Z_U-Z_L).*fu_x_full.*fell_x_full;
  B02_x = -(1.0/(2.0*hbar)^2)*fu_full.*fu_x_full;
  B02_z = (1.0/(2.0*hbar))*Z_L.*fu_x_full.^2;
  
  B1_x = B10_x + B01_x;
  %B1_z = 0
  B2_x = B20_x + B11_x + B02_x;
  B2_z = B20_z + B11_z + B02_z;
  
  C10 = -k2*(2.0/(2.0*hbar))*fell_full;
  C01 = k2*(2.0/(2.0*hbar))*fu_full;
  C20 = k2*(1.0/(2.0*hbar)^2)*fell_full.^2;
  C11 = -k2*(2.0/(2.0*hbar)^2)*fu_full.*fell_full;
  C02 = k2*(1.0/(2.0*hbar)^2)*fu_full.^2;
  
  C1 = C10 + C01;
  C2 = C20 + C11 + C02;
  
  if(n>=1)
    u_x = dxp(un(:,:,n-1+1),alphap,eep,eem,Nx,Nz);
    temp = A1_xx.*u_x;
    Fn = Fn - dxp(temp,alphap,eep,eem,Nx,Nz);
    temp = A1_zx.*u_x;
    Fn = Fn - dz(temp,Dz,2*hbar,Nx,Nz);
    temp = B1_x.*u_x;
    Fn = Fn - temp;
    
    u_z = dz(un(:,:,n-1+1),Dz,2*hbar,Nx,Nz);
    temp = A1_xz.*u_z;
    Fn = Fn - dxp(temp,alphap,eep,eem,Nx,Nz);
    %A1_zz = 0
    %B1_z = 0
    
    temp = C1.*un(:,:,n-1+1);
    Fn = Fn - temp;
    
    Hun = Hun ...
        + (1.0/(2*hbar))*fu.*Vun(:,n-1+1)...
        - (1.0/(2*hbar))*fell.*Vun(:,n-1+1)...
        + (1i*eta/(2*hbar))*fu.*un(:,ell_top,n-1+1)...
        - (1i*eta/(2*hbar))*fell.*un(:,ell_top,n-1+1)...
        + fu_x.*u_x(:,ell_top);
    
    Helln = Helln ...
        + (1.0/(2*hbar))*fu.*Velln(:,n-1+1)...
        - (1.0/(2*hbar))*fell.*Velln(:,n-1+1)...
        + (1i*eta/(2*hbar))*fu.*un(:,ell_bottom,n-1+1)...
        - (1i*eta/(2*hbar))*fell.*un(:,ell_bottom,n-1+1)...
        - fell_x.*u_x(:,ell_bottom);
  end
  
  if(n>=2)
    u_x = dxp(un(:,:,n-2+1),alphap,eep,eem,Nx,Nz);
    temp = A2_xx.*u_x;
    Fn = Fn - dxp(temp,alphap,eep,eem,Nx,Nz);
    temp = A2_zx.*u_x;
    Fn = Fn - dz(temp,Dz,2*hbar,Nx,Nz);
    temp = B2_x.*u_x;
    Fn = Fn - temp;
    
    u_z = dz(un(:,:,n-2+1),Dz,2*hbar,Nx,Nz);
    temp = A2_xz.*u_z;
    Fn = Fn - dxp(temp,alphap,eep,eem,Nx,Nz);
    temp = A2_zz.*u_z;
    Fn = Fn - dz(temp,Dz,2*hbar,Nx,Nz);
    temp = B2_z.*u_z;
    Fn = Fn - temp;
    
    temp = C2.*un(:,:,n-2+1);
    Fn = Fn - temp;
    
    Hun = Hun ...
        + (1.0/(2*hbar))*fu.*fu_x.*u_x(:,ell_top)...
        - (1.0/(2*hbar))*fell.*fu_x.*u_x(:,ell_top)...
        - fu_x.*fu_x.*u_z(:,ell_top);

    Helln = Helln ...
        - (1.0/(2*hbar))*fu.*fell_x.*u_x(:,ell_bottom)...
        + (1.0/(2*hbar))*fell.*fell_x.*u_x(:,ell_bottom)...
        + fell_x.*fell_x.*u_z(:,ell_bottom);
  end
  
  % Solve elliptic equation
  
  Fnhat = zeros(Nx,Nz+1);
  for ell=0:Nz
    Fnhat(:,ell+1) = fft(eem.*Fn(:,ell+1));
  end
  Hunhat = fft(eem.*Hun);
  Hellnhat = fft(eem.*Helln);
  
  for j=1:Nx
    Fnhat_p = Fnhat(j,:).';
    alphaalpha = 1.0;
    betabeta = 0.0;
    gammagamma = k2 - (alphap(j))^2;
    d_a = -1i*eta;
    n_a = -1.0;
    r_a = Hellnhat(j);
    d_b = -1i*eta;
    n_b = 1.0;
    r_b = Hunhat(j);
    uhat_p = solvebvp_colloc(Fnhat_p,alphaalpha,betabeta,gammagamma,...
        (2.0/(z_max-z_min))*Dz,d_a,n_a,r_a,d_b,n_b,r_b);
    temphat(j,:) = uhat_p.';
  end
  
  for ell=0:Nz
    un(:,ell+1,n+1) = eep.*ifft(temphat(:,ell+1));
  end

end

return;
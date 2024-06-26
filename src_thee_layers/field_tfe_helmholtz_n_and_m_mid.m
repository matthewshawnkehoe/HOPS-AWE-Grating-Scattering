function [vnm] = field_tfe_helmholtz_n_and_m_mid(Qunm,Qellnm,hbar,eta,fu,fell,...
    p,alpha,alphap,gamma,gammap,eep,eem,Dz,a,Nx,Nz,N,M,identy)
% field_tfe_helmholtz_m_and_n_mid.m: Computes the approximate solution 
% (vnm) in the middle field.
%
%  Inputs:
%   Qunm: a tensor representing the solution in the upper half of the middle
%        field
%   Qellnm: a tensor representing the solution in the lower half of the middle
%        field
%   hbar: a numerical onstant representing the height above the interface
%         (?)
%   eta: a numerical sontant
%   fu: a test function representing the grating surface in the upper half
%       of the middle field
%   fell: a test function representing the grating surface in the lower half
%       of the middle field
%   p: an integer where tilde_p = (2*pi/d)*p and d is the periodicity of the grating interface
%   alphap: a numerical constant for all wave numbers p
%   gammap: a numerical constant for all wave numbers p
%   eep: phase correction (positive part)
%   eem: phase correction (negative part)
%   Dz: the partial derivative with respect to the z component
%   a: the artificial boundary imposed at the top of the upper layer
%   Nx: the number of discretization points
%   Nz: the number of collocation points
%   N: the maximum number of Taylor orders for the interfacial perturbation
%   M: the maximum number of Taylor orders for the frequency perturbation
%   identy: the identity matrix
%
%  Output:
%   vnm: a tensor representing the approximate solution in the middle field

vnm = zeros(Nx,Nz+1,M+1,N+1);
Qunmhat = zeros(Nx,M+1,N+1);
Qellnmhat = zeros(Nx,M+1,N+1);

k2 = alphap(0+1)^2 + gammap(0+1)^2;

ell_top = 0 + 1;
ell_bottom = Nz + 1;
for n=0:N
  for m=0:M % Needed?
    Qunmhat(:,m+1,n+1) = fft(Qunm(:,m+1,n+1));
    Qellnmhat(:,m+1,n+1) = fft(Qellnm(:,m+1,n+1));
  end
end
fu_x = real(ifft( (1i*p).*fft(fu) ));
fell_x = real(ifft( (1i*p).*fft(fell) ));

ll = [0:Nz]';
z_min = -hbar; z_max = hbar;
% Parameters for solvebvp_colloc
D = (2/(z_max-z_min))*Dz;
D2 = D*D;
D_start = D(1,:);
D_end = D(end,:);
z = ((z_max-z_min)/2.0)*(cos(pi*ll/Nz) - 1.0) + z_max;

fu_full = repmat(fu,1,Nz+1);
fell_full = repmat(fell,1,Nz+1);
fu_x_full = repmat(fu_x,1,Nz+1);
fell_x_full = repmat(fell_x,1,Nz+1);
Tv = T_dno(alpha,alphap,gamma,gammap,k2,Nx,M);
temphat = zeros(Nx,Nz+1);

Z_L = zeros(Nx,Nz+1);
Z_U = zeros(Nx,Nz+1);
for ell=0:Nz
  for j=1:Nx
    Z_L(j,ell+1) = (z(ell+1) - z_min)/(2.0*hbar);
    Z_U(j,ell+1) = (z_max - z(ell+1))/(2.0*hbar);
  end
end

% n = 0 and m = 0

sh = sinh(1i*gammap*hbar)./(1i*gammap);
ch = cosh(1i*gammap*hbar);
aa = -(gammap.^2).*sh - 1i*eta*ch;
bb = ch - 1i*eta*sh;
Bp = (Qunmhat(:,0+1,0+1) + Qellnmhat(:,0+1,0+1))./(2.0*aa);
Cp = (Qunmhat(:,0+1,0+1) - Qellnmhat(:,0+1,0+1))./(2.0*bb);
for ell=0:Nz
  vnm(:,ell+1,0+1,0+1) = ifft( Bp.*cosh(1i*gammap*z(ell+1))...
      + Cp.*sinh(1i*gammap*z(ell+1))./(1i*gammap) );
end

% n > 0 or m > 0

% Preallocate Fnm, Jnm, Hunm, and Hellnm
Fnm = zeros(Nx, Nz+1);
Jnm = zeros(Nx, 1);
%Hunm =  zeros(Nx,M+1,N+1);
Hellnm =  zeros(Nx,M+1,N+1);

for n=0:N
  for m=0:M

    % Form Fnm, Jnm, Hunm, Hellnm
    Fnm(:) = 0;
    Jnm(:) = 0;
    % TODO - Correct initialization of Hunm, Hunm?
    Hunm(:) = 0;
    Hellnm(:) = 0;
    
    % TODO - Verify A, B, and C
  
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
      u_x = dxp(vnm(:,:,m+1,n-1+1),alphap,eep,eem,Nx,Nz); % Correct m indexing?
      temp = A1_xx.*u_x;
      Fnm = Fnm - dxp(temp,alphap,eep,eem,Nx,Nz);
      temp = A1_zx.*u_x;
      Fnm = Fnm - dz(temp,Dz,2*hbar);
      temp = B1_x.*u_x;
      Fnm = Fnm - temp;
    
      u_z = dz(vnm(:,:,m+1,n-1+1),Dz,2*hbar); % Correct m indexing?
      temp = A1_xz.*u_z;
      Fnm = Fnm - dxp(temp,alphap,eep,eem,Nx,Nz);
      %A1_zz = 0
      %B1_z = 0
    
      temp = C1.*vnm(:,:,m+1,n-1+1); % Correct m indexing?
      Fnm = Fnm - temp;

      % Copied from upper field
      temp = 2*1i*alpha.*C1.*u_x;
      Fnm = Fnm - temp;
	  temp = gamma^2.*C1.*vnm(:,:,m+1,n-1+1);
	  Fnm = Fnm - temp;

      % TODO - Fix these
      Hunm = Hunm ...
          - (1.0/(hbar))*fu.*Qunm(:,m+1,n-1+1);
          % - fu_x*u_x(:,ell_top,m+1,n-1+1);
    
      Hellnm = Hellnm ...
          - (1.0/(hbar))*fu.*Qellnm(:,m+1,n-1+1);
    end

    if(m>=1)
 	  u_x = dx(vnm(:,:,m-1+1,n+1),p);
 	  temp = 2*1i*alpha.*u_x;
 	  Fnm  = Fnm - temp;
 	  temp = 2*gamma^2.*vnm(:,:,m-1+1,n+1);
 	  Fnm  = Fnm - temp;

      % TODO - Fix these
      Hunm = Hunm ...
          + (1.0/(hbar))*fell.*Qunm(:,m-1+1,n+1);

      Hellnm = Hellnm ...
          + (1.0/(hbar))*fell.*Qellnm(:,m-1+1,n+1) ... 
          + fell_x.*u_x(:,ell_bottom);
    end

    if(n>=1 && m>=1)
 	  u_x = dx(vnm(:,:,m-1+1,n-1+1),p);
 	  temp = 2*1i*alpha.*C1.*u_x;
 	  Fnm  = Fnm - temp;
 	  temp = 2*gamma^2.*C1.*vnm(:,:,m-1+1,n-1+1);
 	  Fnm  = Fnm - temp;

      % TODO - Fix these
      Hunm = Hunm ...
          + (1.0/(hbar))*fell.*fu_x.*u_x(:,ell_top);
      Hellnm = Hellnm ...
          + (1.0/(hbar))*fu.*fell_x.*u_x(:,ell_bottom);
 	end
  
    if(n>=2)
      u_x = dxp(vnm(:,:,m+1,n-2+1),alphap,eep,eem,Nx,Nz); % Correct m indexing?
      temp = A2_xx.*u_x;
      Fnm = Fnm - dxp(temp,alphap,eep,eem,Nx,Nz);
      temp = A2_zx.*u_x;
      Fnm = Fnm - dz(temp,Dz,2*hbar);
      temp = B2_x.*u_x;
      Fnm = Fnm - temp;
    
      u_z = dz(vnm(:,:,m+1,n-2+1),Dz,2*hbar); % Correct m indexing?
      temp = A2_xz.*u_z;
      Fnm = Fnm - dxp(temp,alphap,eep,eem,Nx,Nz);
      temp = A2_zz.*u_z;
      Fnm = Fnm - dz(temp,Dz,2*hbar);
      temp = B2_z.*u_z;
      Fnm = Fnm - temp;
    
      temp = C2.*vnm(:,:,m+1,n-2+1); % Correct m indexing?
      Fnm = Fnm - temp;
      
      % Copied from upper field
      temp = 2*1i*alpha.*C2.*u_x;
      Fnm = Fnm - temp;
	  temp = gamma^2.*C2.*vnm(:,:,m+1,n-2+1);
	  Fnm = Fnm - temp;
    
      % TODO - Fix these
      Hunm = Hunm ...
          - (1.0/(hbar))*fu.*fu_x.*u_x(:,ell_top)...
          + fu_x.*fu_x.*u_z(:,ell_top);
    end

    if(m>=2)
      u_x = dxp(vnm(:,:,m-2+1,n+1),alphap,eep,eem,Nx,Nz); % Correct m indexing?
 	  temp = gamma^2.*vnm(:,:,m-2+1,n+1);
 	  Fnm = Fnm - temp;

      u_z = dz(vnm(:,:,m-2+1,n+1),Dz,2*hbar); % Correct m indexing?

      Hellnm = Hellnm ...
          - (1.0/(hbar))*fell.*fell_x.*u_x(:,ell_bottom)...
          + fell_x.*fell_x.*u_z(:,ell_bottom);
    end

    if(n>=1 && m>=2)
 	  temp = gamma^2.*C1.*vnm(:,:,m-2+1,n-1+1);
 	  Fnm = Fnm - temp;
    end
 	
    if(n>=2 && m>=1)
 	  u_x = dx(vnm(:,:,m-1+1,n-2+1),p);
 	  temp = 2*1i*alpha.*C2.*u_x;
 	  Fnm = Fnm - temp;
 	  temp = 2*gamma^2.*C2.*vnm(:,:,m-1+1,n-2+1);
 	  Fnm = Fnm - temp;
    end
 	
    if(n>=2 && m>=2)
 	  temp = gamma^2.*C2.*vnm(:,:,m-2+1,n-2+1);
 	  Fnm = Fnm - temp;
    end

    for r=0:m-1
      Jnm = Jnm + ifft( (Tv(:,m-r+1)).*fft(vnm(:,ell_top,r+1,n+1)));
    end
    if(n>=1)
      for r=0:m
        Snm = ifft( (Tv(:,m-r+1)).*fft(vnm(:,ell_top,r+1,n-1+1)) );
        Jnm = Jnm - (1.0/a)*fu.*Snm; %f_u or f_ell?
      end
    end
  
    % Solve elliptic equation
  
    Fnmhat = fft(Fnm);
    Jnmhat = fft(Jnm);
    Hunmhat = fft(Hunm);
    Hellnmhat = fft(Hellnm);
  
    b = Fnmhat.';
    alphaalpha = 1.0;
    betabeta = 0.0;
    gammagamma = k2 - alphap.^2;
    d_min = -1i*eta;
    n_min = -1.0;
    r_min = Hellnmhat(:,m+1,n+1);
    d_max = -1i*eta;
    n_max = 1.0;
    r_max = Hunmhat;

    Vhat = solvebvp_colloc_fast(b,alphaalpha,betabeta,gammagamma,...
             d_min,n_min,r_min,d_max,n_max,r_max,Nx,identy,D,D2,D_start,D_end);  % (2.0/(z_max-z_min))*Dz after gammagamma?

    if((n>0)||(m>0))
      % temphat = Vhat.';
      vnm(:,:,m+1,n+1) = ifft(Vhat); % Correct m indexing?
    end 
  
    % for ell=0:Nz
    %   vnm(:,ell+1,m+1,n+1) = ifft(temphat(:,ell+1)); % Correct m indexing?
    % end

  end
end

return;
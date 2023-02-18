function [wmn] = field_tfe_helmholtz_m_and_n_lf(...
    xi_lf_n_m,f,p,gammapw,alpha,gammaw,Dz,b,Nx,Nz,N,M,identy,alphap)

% MSK 7/30/21: Vectorized (most of) j loop

wmn = zeros(Nx,Nz+1,M+1,N+1);

%k2 = p(0+1)^2 + gammapw(0+1)^2;
k2 = alphap(0+1)^2 + gammapw(0+1)^2;

ell_bottom = Nz + 1;
xi_n_m_hat = fft(xi_lf_n_m);
f_x = real(ifft( (1i*p).*fft(f) ));

ll = [0:Nz]';
z_min = -b; z_max = 0;
% MSK 7/22/21 - Parameters for solvebvp_colloc
D = (2/(z_max-z_min))*Dz;
D2 = D*D;
D_start = D(1,:);
D_end = D(end,:);
% End MSK 7/22/21
tilde_z = cos(pi*ll/Nz);
z = ((z_max-z_min)/2.0)*(tilde_z - 1.0) + z_max;

f_full = repmat(f,1,Nz+1);
f_x_full = repmat(f_x,1,Nz+1);
b_plus_z_full = repmat(b + z.',Nx,1);

What = zeros(Nx,Nz+1);

Tw = T_dno(alpha,p,gammaw,gammapw,k2,Nx,M);

% n=0 and m=0

for ell=0:Nz
  wmn(:,ell+1,0+1,0+1) = ifft( exp(-1i*gammapw*z(ell+1)).*xi_n_m_hat(:,0+1,0+1) );
end

A1_xx = (2.0/b)*f_full;
A1_xz = -(1.0/b)*(b_plus_z_full).*f_x_full;
A1_zx = A1_xz;
%A1_zz = 0;
  
A2_xx = (1.0/b^2)*f_full.^2;
A2_xz = -(1.0/b^2)*(b_plus_z_full).*(f_full.*f_x_full);
A2_zx = A2_xz;
A2_zz = (1.0/b^2)*((b_plus_z_full).^2).*(f_x_full.^2);
  
B1_x = -(1.0/b)*f_x_full;
%B1_z = 0;
  
B2_x = -(1.0/b^2)*f_full.*f_x_full;
B2_z = (1.0/b^2).*(b_plus_z_full).*(f_x_full.^2);
  
S1 = (2.0/b)*f_full;
S2 = (1.0/b^2)*f_full.^2;

for n=0:N
  for m=0:M
    
    % Form Fnm, Qnm
    Fnm = zeros(Nx,Nz+1);
    Qnm = zeros(Nx,1);
    
    if(n>=1)
      w_x = dx(wmn(:,:,m+1,n-1+1),p);
      temp = A1_xx.*w_x;
      Fnm = Fnm - dx(temp,p);
      temp = A1_zx.*w_x;
      Fnm = Fnm - dz(temp,Dz,b);
      temp = B1_x.*w_x;
      Fnm = Fnm - temp;
    
      w_z = dz(wmn(:,:,m+1,n-1+1),Dz,b);
      temp = A1_xz.*w_z;
      Fnm = Fnm - dx(temp,p);
      %A1_zz = 0
      %B1_z = 0
    
      temp = 2*1i*alpha.*S1.*w_x;
      Fnm = Fnm - temp;
	  temp = gammaw^2.*S1.*wmn(:,:,m+1,n-1+1);
	  Fnm = Fnm - temp;
    end
	
 	if(m>=1)
 	  w_x = dx(wmn(:,:,m-1+1,n+1),p);
 	  temp = 2*1i*alpha.*w_x;
 	  Fnm  = Fnm - temp;
 	  temp = 2*gammaw^2.*wmn(:,:,m-1+1,n+1);
 	  Fnm  = Fnm - temp;
 	end
	
 	if(n>=1 && m>=1)
 	  w_x = dx(wmn(:,:,m-1+1,n-1+1),p);
 	  temp = 2*1i*alpha.*S1.*w_x;
 	  Fnm  = Fnm - temp;
 	  temp = 2*gammaw^2.*S1.*wmn(:,:,m-1+1,n-1+1);
 	  Fnm  = Fnm - temp;
 	end
	 
    if(n>=2)
      w_x = dx(wmn(:,:,m+1,n-2+1),p);
      temp = A2_xx.*w_x;
      Fnm = Fnm - dx(temp,p);
      temp = A2_zx.*w_x;
      Fnm = Fnm - dz(temp,Dz,b);
      temp = B2_x.*w_x;
      Fnm = Fnm - temp;
    
      w_z = dz(wmn(:,:,m+1,n-2+1),Dz,b);
      temp = A2_xz.*w_z;
      Fnm = Fnm - dx(temp,p);
      temp = A2_zz.*w_z;
      Fnm = Fnm - dz(temp,Dz,b);
      temp = B2_z.*w_z;
      Fnm = Fnm - temp;
    
      temp = 2*1i*alpha.*S2.*w_x;
      Fnm = Fnm - temp;
	  temp = gammaw^2.*S2.*wmn(:,:,m+1,n-2+1);
	  Fnm = Fnm - temp;
    end
	
    if(m>=2)
 	  temp = gammaw^2.*wmn(:,:,m-2+1,n+1);
 	  Fnm = Fnm - temp;
    end
 	
    if(n>=1 && m>=2)
 	 temp = gammaw^2.*S1.*wmn(:,:,m-2+1,n-1+1);
 	 Fnm = Fnm - temp;
    end
 	
    if(n>=2 && m>=1)
 	 w_x = dx(wmn(:,:,m-1+1,n-2+1),p);
 	 temp = 2*1i*alpha.*S2.*w_x;
 	 Fnm = Fnm - temp;
 	 temp = 2*gammaw^2.*S2.*wmn(:,:,m-1+1,n-2+1);
 	 Fnm = Fnm - temp;
    end
 	
    if(n>=2 && m>=2)
 	 temp = gammaw^2.*S2.*wmn(:,:,m-2+1,n-2+1);
 	 Fnm = Fnm - temp;
    end

    for r=0:m-1
      Qnm = Qnm - ifft( (Tw(:,m-r+1)).*fft(wmn(:,ell_bottom,r+1,n+1)));
    end
    if(n>=1)
      for r=0:m
        Snm = ifft( (Tw(:,m-r+1)).*fft(wmn(:,ell_bottom,r+1,n-1+1)) );
        Qnm = Qnm - (1.0/b)*f.*Snm;
      end
    end
    
    % Solve elliptic equation
    
    Fnmhat = fft(Fnm);
    Qnmhat = fft(Qnm);
    
    b_lower = Fnmhat.';
    alphaalpha = 1.0;
    betabeta = 0.0;
    %gammagamma = gammaw*gammaw - p.^2 - 2*alpha.*p;
    %gammagamma = gammapw.^2;
    gammagamma = k2 - alphap.^2;
    d_min = 1i*gammapw;
    n_min = 1.0;
    r_min = Qnmhat;
    d_max = 1.0;
    n_max = 0.0;
    r_max = xi_n_m_hat(:,m+1,n+1);
    
    What = solvebvp_colloc_fast_lf(b_lower,alphaalpha,betabeta, ...
        gammagamma,d_min,n_min,r_min,d_max,n_max,r_max,Nx,identy,D,D2,...
        D_start,D_end);
    
    if((n>0)||(m>0))
      wmn(:,:,m+1,n+1) = ifft(What);
    end
  end
end

return;

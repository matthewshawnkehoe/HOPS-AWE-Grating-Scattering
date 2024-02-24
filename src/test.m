clear all; close all; clc;

% Constants
a = 1.0;
alpha = 0;
gamma = 1.5;
Nx = 32;
Nz = 32;
N = 20;
M = 20;

% N-Dimensional Test Matrices - these have a lot more meaning in a larger
% program where they are being called
p = rand(Nx,1);
real_part = rand(Nx, Nz+1, N+1, M+1);
imaginary_part = rand(Nx, Nz+1, N+1, M+1);
unm = complex(real_part, imaginary_part);

A1_xx = rand(Nx, Nz+1);
A1_xz = rand(Nx, Nz+1);
A1_zx = rand(Nx, Nz+1);
  
A2_xx = rand(Nx, Nz+1);
A2_xz = rand(Nx, Nz+1);
A2_zx = rand(Nx, Nz+1);
A2_zz = rand(Nx, Nz+1);
  
B1_x = rand(Nx, Nz+1);
  
B2_x = rand(Nx, Nz+1);
B2_z = rand(Nx, Nz+1);

Dz = rand(Nx+1, Nz+1);

S1 = rand(Nx, Nz+1);
S2 = rand(Nx, Nz+1);

% Preallocate Fnm
Fnm = zeros(Nx, Nz+1);

% Preallocate u_x and u_z
u_x = zeros(Nx, Nz+1, M+1, N+1);
u_z = zeros(Nx, Nz+1, M+1, N+1);

% Double for loop
for n=0:N
  for m=0:M
    
    Fnm(:) = 0;
    
    if(n>=1)
      u_x = dx(unm(:,:,m+1,n-1+1),p); % Expensive to call repeatedly in nested loop
      temp = A1_xx.*u_x;
      Fnm = Fnm - dx(temp,p);
      temp = A1_zx.*u_x;
      Fnm = Fnm - dz(temp,Dz,a);
      temp = B1_x.*u_x;
      Fnm = Fnm - temp;
    
      u_z = dz(unm(:,:,m+1,n-1+1),Dz,a); % Expensive to call repeatedly in nested loop
      temp = A1_xz.*u_z;
      Fnm = Fnm - dx(temp,p);
      temp = 2*1i*alpha.*S1.*u_x;
      Fnm = Fnm - temp;
	  temp = gamma^2.*S1.*unm(:,:,m+1,n-1+1);
	  Fnm = Fnm - temp;
    end
	
 	if(m>=1)
 	  u_x = dx(unm(:,:,m-1+1,n+1),p);
 	  temp = 2*1i*alpha.*u_x;
 	  Fnm  = Fnm - temp;
 	  temp = 2*gamma^2.*unm(:,:,m-1+1,n+1);
 	  Fnm  = Fnm - temp;
 	end
	
 	if(n>=1 && m>=1)
 	  u_x = dx(unm(:,:,m-1+1,n-1+1),p);
 	  temp = 2*1i*alpha.*S1.*u_x;
 	  Fnm  = Fnm - temp;
 	  temp = 2*gamma^2.*S1.*unm(:,:,m-1+1,n-1+1);
 	  Fnm  = Fnm - temp;
 	end
	 
    if(n>=2)
      u_x = dx(unm(:,:,m+1,n-2+1),p);
      temp = A2_xx.*u_x;
      Fnm = Fnm - dx(temp,p);
      temp = A2_zx.*u_x;
      Fnm = Fnm - dz(temp,Dz,a);
      temp = B2_x.*u_x;
      Fnm = Fnm - temp;
    
      u_z = dz(unm(:,:,m+1,n-2+1),Dz,a);
      temp = A2_xz.*u_z;
      Fnm = Fnm - dx(temp,p);
      temp = A2_zz.*u_z;
      Fnm = Fnm - dz(temp,Dz,a);
      temp = B2_z.*u_z;
      Fnm = Fnm - temp;
    
      temp = 2*1i*alpha.*S2.*u_x;
      Fnm = Fnm - temp;
	  temp = gamma^2.*S2.*unm(:,:,m+1,n-2+1);
	  Fnm = Fnm - temp;
    end
	
    if(m>=2)
 	  temp = gamma^2.*unm(:,:,m-2+1,n+1);
 	  Fnm = Fnm - temp;
    end
 	
    if(n>=1 && m>=2)
 	 temp = gamma^2.*S1.*unm(:,:,m-2+1,n-1+1);
 	 Fnm = Fnm - temp;
    end
 	
    if(n>=2 && m>=1)
 	 u_x = dx(unm(:,:,m-1+1,n-2+1),p);
 	 temp = 2*1i*alpha.*S2.*u_x;
 	 Fnm = Fnm - temp;
 	 temp = 2*gamma^2.*S2.*unm(:,:,m-1+1,n-2+1);
 	 Fnm = Fnm - temp;
    end
 	
    if(n>=2 && m>=2)
 	 temp = gamma^2.*S2.*unm(:,:,m-2+1,n-2+1);
 	 Fnm = Fnm - temp;
    end
  end
end

display(Fnm)

return;

function [u_z] = dz(u,Dz,b)

u_z = ((2.0/b)*Dz*u.').';

return;

function [u_x] = dx(u,p)

u_x = ifft((1i*p).*fft(u));

return;
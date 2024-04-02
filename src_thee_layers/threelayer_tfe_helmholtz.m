function [U_nm,Vu_nm,Vell_nm,W_nm] ...
  = threelayer_tfe_helmholtz(zeta_n,psi_n,theta_n,mu_n,...
  hbar,eta,fu,fell,tau2,sigma2,p,alphap,gamma_up,gamma_vp,gamma_wp,...
  eep,eem,a,b,Nx,Nz,N)

% TODO - Adjust to your problem

[Dz,z] = cheb(Nz);

% hbar is the "half-height" of the layer

U_nm = zeros(Nx,N+1);
Vu_nm = zeros(Nx,N+1);
Vell_nm = zeros(Nx,N+1);
W_nm = zeros(Nx,N+1);
Q_U_rs = zeros(Nx,N+1,N+1);
Ru_V_rs = zeros(Nx,N+1,N+1);
Rell_V_rs = zeros(Nx,N+1,N+1);
S_W_rs = zeros(Nx,N+1,N+1);

Q0 = (-1i*gamma_up+1i*eta)./(-1i*gamma_up-1i*eta);
S0 = (-1i*gamma_wp+1i*eta)./(-1i*gamma_wp-1i*eta);
sh = sinh(1i*gamma_vp*hbar)./(1i*gamma_vp);
ch = cosh(1i*gamma_vp*hbar);
aa = -(gamma_vp.^2).*sh - 1i*eta*ch;
bb = ch - 1i*eta*sh;
Ruu0 = (1.0/2.0)*( conj(aa)./aa + conj(bb)./bb );
Ruell0 = (1.0/2.0)*( conj(aa)./aa - conj(bb)./bb );
Rellu0 = Ruell0;
Rellell0 = Ruu0;

for n=0:N
  zetahat = fft(eem.*zeta_n(:,n+1));
  psihat = fft(eem.*psi_n(:,n+1));
  thetahat = fft(eem.*theta_n(:,n+1));
  muhat = fft(eem.*mu_n(:,n+1));
  Y_Uhat = -2*1i*eta*zetahat;
  Y_Vuhat = -2*psihat;
  Y_Vellhat = -2*1i*eta*thetahat;
  Y_What = -2*muhat;

  for m=0:n-1
    YY = -Q_U_rs(:,n-m+1,m+1) + Ru_V_rs(:,n-m+1,m+1);
    Y_Uhat = Y_Uhat - fft(eem.*YY);
    YY = Q_U_rs(:,n-m+1,m+1) + tau2*Ru_V_rs(:,n-m+1,m+1);
    Y_Vuhat = Y_Vuhat - fft(eem.*YY);
    YY = -Rell_V_rs(:,n-m+1,m+1) + S_W_rs(:,n-m+1,m+1);
    Y_Vellhat = Y_Vellhat - fft(eem.*YY);
    YY = Rell_V_rs(:,n-m+1,m+1) + sigma2*S_W_rs(:,n-m+1,m+1);
    Y_What = Y_What - fft(eem.*YY);
  end

  Uhat = zeros(Nx,1);
  Vuhat = zeros(Nx,1);
  Vellhat = zeros(Nx,1);
  What = zeros(Nx,1);

  for j=1:Nx
    MM = [1.0-Q0(j) -1.0+Ruu0(j) Ruell0(j) 0.0;...
        1.0+Q0(j) tau2*(1.0+Ruu0(j)) tau2*Ruell0(j) 0.0;...
        0.0 -Rellu0(j) 1.0-Rellell0(j) -1.0+S0(j);...
        0.0 Rellu0(j) 1.0+Rellell0(j) sigma2*(1.0+S0(j))];
    bb = [Y_Uhat(j);Y_Vuhat(j);Y_Vellhat(j);Y_What(j)];
    xx = MM\bb;
    Uhat(j) = xx(1);
    Vuhat(j) = xx(2);
    Vellhat(j) = xx(3);
    What(j) = xx(4);
  end

  U_nm(:,n+1) = eep.*ifft(Uhat);
  Vu_nm(:,n+1) = eep.*ifft(Vuhat);
  Vell_nm(:,n+1) = eep.*ifft(Vellhat);
  W_nm(:,n+1) = eep.*ifft(What);

  % Compute and store Q_r[U_s]
  s = n;
  xi = zeros(Nx,N-s+1);
  xi(:,0+1) = U_nm(:,s+1);
  unm = field_tfe_helmholtz_m_and_n(xi,eta,fu,...
      p,alphap,gamma_up,eep,eem,Dz,a,Nx,Nz,N-s);
  Qnm = iio_tfe_helmholtz_upper(unm,eta,fu,...
      p,alphap,gamma_up,eep,eem,Dz,a,Nx,Nz,N-s);
  for r=0:N-s
    Q_U_rs(:,r+1,s+1) = Qnm(:,r+1);
  end

  % Compute and store Ru_r[V_s], Rell_r[V_s]
  s = n;
  xi = zeros(Nx,N-s+1);
  zeta = zeros(Nx,N-s+1);
  xi(:,0+1) = Vu_nm(:,s+1);
  zeta(:,0+1) = Vell_nm(:,s+1);
  vnm = field_tfe_helmholtz_m_and_n_mid(xi,zeta,hbar,eta,fu,fell,...
      p,alphap,gamma_vp,eep,eem,Dz,a,Nx,Nz,N-s);
  [Runm,Rellnm] = iio_tfe_helmholtz_middle(vnm,hbar,eta,fu,fell,...
      p,alphap,gamma_vp,eep,eem,Dz,a,Nx,Nz,N-s);
  for r=0:N-s
    Ru_V_rs(:,r+1,s+1) = Runm(:,r+1);
    Rell_V_rs(:,r+1,s+1) = Rellnm(:,r+1);
  end

  % Compute and store S_r[W_s]
  s = n;
  xi = zeros(Nx,N-s+1);
  xi(:,0+1) = W_nm(:,s+1);
  wnm = field_tfe_helmholtz_m_and_n_lf(xi,eta,fell,...
      p,alphap,gamma_wp,eep,eem,Dz,b,Nx,Nz,N-s);
  Snm = iio_tfe_helmholtz_lower(wnm,eta,fell,...
      p,alphap,gamma_wp,eep,eem,Dz,b,Nx,Nz,N-s);
  for r=0:N-s
    S_W_rs(:,r+1,s+1) = Snm(:,r+1);
  end
  
end

return;
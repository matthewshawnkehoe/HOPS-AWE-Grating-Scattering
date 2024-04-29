function [U_n_m,Vu_n_m,Vell_n_m,W_n_m] ...
  = threelayer_tfe_helmholtz(zeta_n,psi_n,theta_n,mu_n,...
  hbar,eta,fu,fell,tau2,sigma2,p,alpha,gamma_u,gamma_v,gamma_w,...
  alphap,gamma_up,gamma_vp,gamma_wp,eep,eem,a,b,Nx,Nz,N,M,identy)

% TODO - Adjust to your problem

[Dz,z] = cheb(Nz);

% hbar is the "half-height" of the layer

U_n_m = zeros(Nx,N+1,M+1);
Vu_n_m = zeros(Nx,N+1,M+1);
Vell_n_m = zeros(Nx,N+1,M+1);
W_n_m = zeros(Nx,N+1,M+1);
Q_U_rs = zeros(Nx,N+1,M+1);
Ru_V_rs = zeros(Nx,N+1,M+1);
Rell_V_rs = zeros(Nx,N+1,M+1);
S_W_rs = zeros(Nx,N+1,M+1);

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
  for m=0:M
    zetahat = fft(eem.*zeta_n(:,n+1,m+1));
    psihat = fft(eem.*psi_n(:,n+1,m+1));
    thetahat = fft(eem.*theta_n(:,n+1,m+1));
    muhat = fft(eem.*mu_n(:,n+1,m+1));
    Y_Uhat = -2*1i*eta*zetahat;
    Y_Vuhat = -2*psihat;
    Y_Vellhat = -2*1i*eta*thetahat;
    Y_What = -2*muhat;

    for k=0:n-1
      YY = -Q_U_rs(:,n-k+1,k+1) + Ru_V_rs(:,n-k+1,k+1);
      Y_Uhat = Y_Uhat - fft(eem.*YY);
      YY = Q_U_rs(:,n-k+1,k+1) + tau2*Ru_V_rs(:,n-k+1,k+1);
      Y_Vuhat = Y_Vuhat - fft(eem.*YY);
      YY = -Rell_V_rs(:,n-k+1,k+1) + S_W_rs(:,n-k+1,k+1);
      Y_Vellhat = Y_Vellhat - fft(eem.*YY);
      YY = Rell_V_rs(:,n-k+1,k+1) + sigma2*S_W_rs(:,n-k+1,k+1);
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

    U_n_m(:,n+1,m+1) = eep.*ifft(Uhat);
    Vu_n_m(:,n+1,m+1) = eep.*ifft(Vuhat);
    Vell_n_m(:,n+1,m+1) = eep.*ifft(Vellhat);
    W_n_m(:,n+1,m+1) = eep.*ifft(What);

    % Compute and store Q_r[U_s]?
    % Compute and store Q_{p,r}[U_{q,s}]
    q = n; s = m;
    xi_n_m = zeros(Nx,M-s+1,N-q+1);
    xi_n_m(:,0+1,0+1) = U_n_m(:,s+1,q+1);
    %eta_n_m = ...
    u_n_m = field_tfe_helmholtz_m_and_n(xi_n_m,fu,p,gamma_up,...
              alpha,gamma_u,Dz,a,Nx,Nz,N-q,M-s,identy,alphap);
    Q_n_m = iio_tfe_helmholtz_upper(u_n_m,...
              f,p,Dz,a,Nx,Nz,N-q,M-s);
    for p=0:N-q
      for r=0:M-s
        Q_U_pqrs(:,r+1,s+1,p+1,q+1) = Q_n_m(:,r+1,p+1);
      end
    end
    % s = n;
    % xi = zeros(Nx,N-s+1);
    % xi(:,0+1) = U_n_m(:,s+1);
    % unm = field_tfe_helmholtz_m_and_n(xi,eta,fu,...
    %     p,alphap,gamma_up,eep,eem,Dz,a,Nx,Nz,N-s);
    % Qnm = iio_tfe_helmholtz_upper(unm,eta,fu,...
    %     p,alphap,gamma_up,eep,eem,Dz,a,Nx,Nz,N-s);
    % for r=0:N-s
    %   Q_U_rs(:,r+1,s+1) = Qnm(:,r+1);
    % end

    % Compute and store Ru_r[V_s], Rell_r[V_s]?
    % Compute and store Ru_{p,r}[V_{q,s}], Rell_{p,r}[V_{q,s}]
    q = n; s = m;
    xi_n_m = zeros(Nx,M-s+1,N-q+1);
    zeta_n_m = zeros(Nx,M-s+1,N-q+1);
    xi_n_m(:,0+1) = Vu_n_m(:,s+1,q+1);
    zeta_n_m(:,0+1) = Vell_n_m(:,s+1,q+1);
    v_n_m = field_tfe_helmholtz_m_and_n_mid(xi,zeta,hbar,eta,fu,fell,...
                        p,alphap,gamma_vp,eep,eem,Dz,a,Nx,Nz,N-s);
    [Runm,Rellnm] = iio_tfe_helmholtz_middle(vnm,hbar,eta,fu,fell,...
                        p,alphap,gamma_vp,eep,eem,Dz,a,Nx,Nz,N-s);
    for p=0:N-q
      for r=0:M-s
        Ru_V_pqrs(:,r+1,s+1,p+1,q+1) = Runm(:,r+1,p+1);
        Rell_V_pqrs(:,r+1,s+1,p+1,q+1) = Rellnm(:,r+1,p+1);
      end
    end

    % s = n;
    % xi = zeros(Nx,N-s+1);
    % zeta = zeros(Nx,N-s+1);
    % xi(:,0+1) = Vu_n_m(:,s+1);
    % zeta(:,0+1) = Vell_n_m(:,s+1);
    % vnm = field_tfe_helmholtz_m_and_n_mid(xi,zeta,hbar,eta,fu,fell,...
    %     p,alphap,gamma_vp,eep,eem,Dz,a,Nx,Nz,N-s);
    % [Runm,Rellnm] = iio_tfe_helmholtz_middle(vnm,hbar,eta,fu,fell,...
    %     p,alphap,gamma_vp,eep,eem,Dz,a,Nx,Nz,N-s);
    % for r=0:N-s
    %   Ru_V_rs(:,r+1,s+1) = Runm(:,r+1);
    %   Rell_V_rs(:,r+1,s+1) = Rellnm(:,r+1);
    % end

    % Compute and store S_r[W_s]?
    % Compute and store S_{p,r}[W_{q,s}]
    q = n; s = m;
    xi_n_m = zeros(Nx,M-s+1,N-q+1);
    xi_n_m(:,0+1,0+1) = W_n_m(:,s+1,q+1);
    %eta_n_m = ...
    wnm = field_tfe_helmholtz_m_and_n_lf(xi_n_m,eta,fell,...
            p,alphap,gamma_wp,eep,eem,Dz,b,Nx,Nz,N-s);
    Snm = iio_tfe_helmholtz_lower(wnm,eta,fell,...
            p,alphap,gamma_wp,eep,eem,Dz,b,Nx,Nz,N-s);
    for p=0:N-q
      for r=0:M-s
        S_W_pqrs(:,r+1,s+1,p+1,q+1) = Snm(:,r+1,p+1);
      end
    end

    % s = n;
    % xi = zeros(Nx,N-s+1);
    % xi(:,0+1) = W_n_m(:,s+1);
    % wnm = field_tfe_helmholtz_m_and_n_lf(xi,eta,fell,...
    %     p,alphap,gamma_wp,eep,eem,Dz,b,Nx,Nz,N-s);
    % Snm = iio_tfe_helmholtz_lower(wnm,eta,fell,...
    %     p,alphap,gamma_wp,eep,eem,Dz,b,Nx,Nz,N-s);
    % for r=0:N-s
    %   S_W_rs(:,r+1,s+1) = Snm(:,r+1);
    % end

  end
end

return;
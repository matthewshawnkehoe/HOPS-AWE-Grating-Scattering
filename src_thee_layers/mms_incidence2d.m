function [U_nm,Utilde_nm,U,Utilde,V_u_nm,Vtilde_u_nm,V_u,Vtilde_u,...
    V_ell_nm,Vtilde_ell_nm,V_ell,Vtilde_ell,W_nm,Wtilde_nm,W,Wtilde] ...
    = mms_incidence2d(hbar,eta,fu,fu_x,fell,fell_x,x,...
    alphap,gamma_up,gamma_vp,gamma_wp,Nx,N,M,Eps,Delta)

% Upper Layer

f = fu; f_x = fu_x;

A_u_r = -3.0; r = 2; alphap_r = alphap(r+1); gamma_up_r = gamma_up(r+1);
xi_u_nm = zeros(Nx,N+1);
nu_u_nm = zeros(Nx,N+1);
f_nm1 = ones(Nx,1);
f_nm = ones(Nx,1);
xi_u_nm(:,0+1) = A_u_r*exp(1i*alphap_r*x);
nu_u_nm(:,0+1) = (-1i*gamma_up_r)*A_u_r*exp(1i*alphap_r*x);
for n=1:N
  if(n>1)
    f_nm1 = f.*f_nm1/(n-1);
  end
  f_nm = f.*f_nm/n;
  xi_u_nm(:,n+1) = A_u_r*exp(1i*alphap_r*x).*(f_nm*(1i*gamma_up_r)^n);
  nu_u_nm(:,n+1) = (-1i*gamma_up_r)*A_u_r*exp(1i*alphap_r*x).*(f_nm*(1i*gamma_up_r)^n)...
      + (1i*alphap_r)*A_u_r*f_x.*exp(1i*alphap_r*x).*(f_nm1*(1i*gamma_up_r)^(n-1));
end
U_nm = nu_u_nm - 1i*eta*xi_u_nm;
Utilde_nm = nu_u_nm + 1i*eta*xi_u_nm;

xi_u = A_u_r*exp(1i*alphap_r*x).*exp(1i*gamma_up_r*Eps.*f);
nu_u = (-1i*gamma_up_r+1i*alphap_r*Eps.*f_x).*xi_u;
U = nu_u - 1i*eta*xi_u;
Utilde = nu_u + 1i*eta*xi_u;

% Lower Layer

f = fell; f_x = fell_x;

A_w_r = 4.0; r = 3; alphap_r = alphap(r+1); gamma_wp_r = gamma_wp(r+1);
xi_w_nm = zeros(Nx,N+1);
nu_w_nm = zeros(Nx,N+1);
f_nm1 = ones(Nx,1);
f_nm = ones(Nx,1);
xi_w_nm(:,0+1) = A_w_r*exp(1i*alphap_r*x);
nu_w_nm(:,0+1) = (-1i*gamma_wp_r)*A_w_r*exp(1i*alphap_r*x);
for n=1:N
  if(n>1)
    f_nm1 = f.*f_nm1/(n-1);
  end
  f_nm = f.*f_nm/n;
  xi_w_nm(:,n+1) = A_w_r*exp(1i*alphap_r*x).*(f_nm*(-1i*gamma_wp_r)^n);
  nu_w_nm(:,n+1) = (-1i*gamma_wp_r)*A_w_r*exp(1i*alphap_r*x).*(f_nm*(-1i*gamma_wp_r)^n)...
      - (1i*alphap_r)*A_w_r*f_x.*exp(1i*alphap_r*x).*(f_nm1*(-1i*gamma_wp_r)^(n-1));
end
W_nm = nu_w_nm - 1i*eta*xi_w_nm;
Wtilde_nm = nu_w_nm + 1i*eta*xi_w_nm;

xi_w = A_w_r*exp(1i*alphap_r*x).*exp(-1i*gamma_wp_r*Eps.*f);
nu_w = (-1i*gamma_wp_r-1i*alphap_r*Eps.*f_x).*xi_w;
W = nu_w - 1i*eta*xi_w;
Wtilde = nu_w + 1i*eta*xi_w;

% Middle Layer

B_v_r = -exp(1); C_v_r = pi;
r = 3; alphap_r = alphap(r+1); gamma_vp_r = gamma_vp(r+1);
xi_v_nm = zeros(Nx,N+1);
zeta_v_nm = zeros(Nx,N+1);
nu_v_nm = zeros(Nx,N+1);
psi_v_nm = zeros(Nx,N+1);
fu_nm1 = ones(Nx,1);
fell_nm1 = ones(Nx,1);
fu_nm = ones(Nx,1);
fell_nm = ones(Nx,1);
eigh = exp(1i*gamma_vp_r*hbar);
eigmh = exp(1i*gamma_vp_r*(-hbar));
xi_v_nm(:,0+1) = (B_v_r*eigh+C_v_r*eigmh)*exp(1i*alphap_r*x);
zeta_v_nm(:,0+1) = (B_v_r*eigmh+C_v_r*eigh)*exp(1i*alphap_r*x);
nu_v_nm(:,0+1) = (1i*gamma_vp_r*B_v_r*eigh-1i*gamma_vp_r*C_v_r*eigmh)...
    *exp(1i*alphap_r*x);
psi_v_nm(:,0+1) = (-1i*gamma_vp_r*B_v_r*eigmh+1i*gamma_vp_r*C_v_r*eigh)...
    *exp(1i*alphap_r*x);
for n=1:N
  if(n>1)
    fu_nm1 = fu.*fu_nm1/(n-1);
    fell_nm1 = fell.*fell_nm1/(n-1);
  end
  fu_nm = fu.*fu_nm/n;
  fell_nm = fell.*fell_nm/n;
  
  temp1 = B_v_r*eigh*exp(1i*alphap_r*x).*(fu_nm*(1i*gamma_vp_r)^n);
  temp2 = C_v_r*eigmh*exp(1i*alphap_r*x).*(fu_nm*(-1i*gamma_vp_r)^n);
  xi_v_nm(:,n+1) = temp1 + temp2;
  
  temp3 = B_v_r*eigh*exp(1i*alphap_r*x).*(fu_nm1*(1i*gamma_vp_r)^(n-1));
  temp4 = C_v_r*eigmh*exp(1i*alphap_r*x).*(fu_nm1*(-1i*gamma_vp_r)^(n-1));
  nu_v_nm(:,n+1) = (1i*gamma_vp_r)*temp1 + (-1i*gamma_vp_r)*temp2...
      - (1i*alphap_r)*fu_x.*temp3 - (1i*alphap_r)*fu_x.*temp4;
  
  temp5 = B_v_r*eigmh*exp(1i*alphap_r*x).*(fell_nm*(1i*gamma_vp_r)^n);
  temp6 = C_v_r*eigh*exp(1i*alphap_r*x).*(fell_nm*(-1i*gamma_vp_r)^n);
  zeta_v_nm(:,n+1) = temp5 + temp6;
  
  temp7 = B_v_r*eigmh*exp(1i*alphap_r*x).*(fell_nm1*(1i*gamma_vp_r)^(n-1));
  temp8 = C_v_r*eigh*exp(1i*alphap_r*x).*(fell_nm1*(-1i*gamma_vp_r)^(n-1));
  psi_v_nm(:,n+1) = (-1i*gamma_vp_r)*temp5 + (1i*gamma_vp_r)*temp6...
      + (1i*alphap_r)*fell_x.*temp7 + (1i*alphap_r)*fell_x.*temp8;
end

V_u_nm = nu_v_nm - 1i*eta*xi_v_nm;
Vtilde_u_nm = nu_v_nm + 1i*eta*xi_v_nm;
V_ell_nm = psi_v_nm - 1i*eta*zeta_v_nm;
Vtilde_ell_nm = psi_v_nm + 1i*eta*zeta_v_nm;

temp1 = B_v_r*exp(1i*alphap_r*x).*exp(1i*gamma_vp_r*(hbar+Eps.*fu));
temp2 = C_v_r*exp(1i*alphap_r*x).*exp(-1i*gamma_vp_r*(hbar+Eps.*fu));
xi_v = temp1 + temp2;
nu_v = (1i*gamma_vp_r)*temp1 + (-1i*gamma_vp_r)*temp2...
    - (1i*alphap_r)*(Eps.*fu_x).*temp1 - (1i*alphap_r)*(Eps.*fu_x).*temp2;
V_u = nu_v - 1i*eta*xi_v;
Vtilde_u = nu_v + 1i*eta*xi_v;
temp5 = B_v_r*exp(1i*alphap_r*x).*exp(1i*gamma_vp_r*(-hbar+Eps.*fell));
temp6 = C_v_r*exp(1i*alphap_r*x).*exp(-1i*gamma_vp_r*(-hbar+Eps.*fell));
zeta_v = temp5 + temp6;
psi_v = -(1i*gamma_vp_r)*temp5 - (-1i*gamma_vp_r)*temp6...
    + (1i*alphap_r)*(Eps.*fell_x).*temp5 + (1i*alphap_r)*(Eps.*fell_x).*temp6;
V_ell = psi_v - 1i*eta*zeta_v;
Vtilde_ell = psi_v + 1i*eta*zeta_v;

return;
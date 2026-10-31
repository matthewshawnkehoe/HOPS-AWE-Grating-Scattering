function [U_n_m,W_n_m,ubar_n_m,wbar_n_m] ...
    = two_layer_solve_coupled(tau2,zeta_n_m,psi_n_m,...
    gamma_u_bar_p,gamma_w_bar_p,N,Nx,f,f_x,pp,alpha_bar,...
    gamma_u_bar,gamma_w_bar,Dz,a,b,Nz,M,identy,alpha_bar_p)
% two_layer_solve_fast.m: Efficiently solves the two-layer problem in the 
% upper and lower fields through storing and reusing coefficients and 
% dynamic programming.
%
%  Inputs:
%   tau2: a numerical constant representing TE or TM mode
%   zeta_n_m: a tensor representing the Dirichlet data
%   psi_n_m: a tensor representing the Neumannn data
%   gamma_u_bar_p: a numerical parameter calculated at every wave number p in the upper field
%   gamma_w_bar_p: a numerical parameter calculated at every wave number p in the lower field
%   N: the maximum number of Taylor orders for the interfacial perturbation
%   Nx: the number of discretization points 
%   f: a test function representing the grating surface
%   f_x: the derivative of the test function in the x component
%   pp: an integer where tilde_p = (2*pi/d)*p and d is the periodicity of the grating interface
%   alpha_bar: a numerical constant
%   gamma_u_bar: a numerical constant in the upper field 
%   gamma_w_bar: a numerical constant in the lower field 
%   Dz: the partial derivative with respect to the z component
%   a: the artificial boundary imposed at the top of the upper layer
%   b: the artificial boundary imposed at the bottom of the lower layer
%   Nz: the number of collocation points
%   M: the maximum number of Taylor orders for the frequency perturbation
%   identy: the identity matrix
%   alpha_bar_p: a numerical constant at all wave numbers p
%
%  Outputs:
%   U_n_m: Approximate solution at the surface z=g(x) in the upper field
%   W_n_m: Approximate solution at the surface z=g(x) in the lower field
%   ubar_n_m: Approximate solution from the upper field solver at the top
%             collocation point
%   wbar_n_m: Approximate solution from the lower field solver at the lower
%             collocation point
% MSK 10/31/26

U_n_m = zeros(Nx,M+1,N+1);
W_n_m = zeros(Nx,M+1,N+1);
ubar_n_m = zeros(Nx,M+1,N+1);
wbar_n_m = zeros(Nx,M+1,N+1);

Lu = layer_setup(true, f,pp,gamma_u_bar_p,alpha_bar,gamma_u_bar,Dz,a,Nx,Nz,N,M,alpha_bar_p);
Lw = layer_setup(false,f,pp,gamma_w_bar_p,alpha_bar,gamma_w_bar,Dz,b,Nx,Nz,N,M,alpha_bar_p);

for n=0:N
  Lu.V(:,:,:,mod(n,3)+1) = 0;            % only three n-levels of the volume field are kept
  Lw.V(:,:,:,mod(n,3)+1) = 0;
  for m=0:M
    [Lu,Gp] = partial(Lu,n,m);           % DNO of the particular field = sum over lower orders
    [Lw,Jp] = partial(Lw,n,m);
    Q = zeta_n_m(:,m+1,n+1);
    R = -psi_n_m(:,m+1,n+1) - Gp - tau2*Jp;
    [U_n_m(:,m+1,n+1),W_n_m(:,m+1,n+1)] = AInverse(Q,R,gamma_u_bar_p,gamma_w_bar_p,Nx,tau2);
    Lu = complete(Lu,n,m,U_n_m(:,m+1,n+1));
    Lw = complete(Lw,n,m,W_n_m(:,m+1,n+1));
    ubar_n_m(:,m+1,n+1) = Lu.tr_bc(:,m+1,n+1);
    wbar_n_m(:,m+1,n+1) = Lw.tr_bc(:,m+1,n+1);
  end
end
end

% -------------------------------------------------------------------------------------
function S = layer_setup(upper,f,p,gammap,alpha,gamma,Dz,L,Nx,Nz,N,M,alphap)
% Everything that does not depend on the data (coefficients, inverted collocation
% matrices, T_dno), exactly as in field_tfe_helmholtz_m_and_n(_lf).m.
S.upper = upper; S.L = L; S.p = p; S.f = f; S.alpha = alpha; S.g2 = gamma^2;
S.Nx = Nx; S.Nz = Nz;
k2 = alphap(1)^2 + gammap(1)^2;
S.fx_r = real(ifft((1i*p).*fft(f)));      % f_x in the field solver (real)
S.fx_c = ifft((1i*p).*fft(f));            % f_x in the DNO (complex)
ll = (0:Nz)';
if upper, z_min = 0; z_max = L; else, z_min = -L; z_max = 0; end
z = ((z_max-z_min)/2.0)*(cos(pi*ll/Nz) - 1.0) + z_max;
D = (2/(z_max-z_min))*Dz; D2 = D*D;
S.DzLt = ((2.0/L)*Dz).';                  % u*DzLt == dz(u,Dz,L)
F = repmat(f,1,Nz+1); FX = repmat(S.fx_r,1,Nz+1);
if upper
  h = repmat(L - z.',Nx,1);
  S.A1_xx = -(2.0/L)*F;            S.A1_xz = -(1.0/L)*h.*FX;
  S.A2_xx = (1.0/L^2)*F.^2;        S.A2_xz = (1.0/L^2)*h.*(F.*FX);
  S.A2_zz = (1.0/L^2)*(h.^2).*(FX.^2);
  S.B1_x = (1.0/L)*FX;             S.B2_x = -(1.0/L^2)*F.*FX;
  S.B2_z = -(1.0/L^2)*h.*(FX.^2);
  S.S1 = -(2.0/L)*F;               S.S2 = (1.0/L^2)*F.^2;
  S.ell_bc = 1; S.ell_if = Nz+1; S.sgn_T = +1;
  S.phase = exp(1i*gammap*z.');    % (Nx, Nz+1)
else
  h = repmat(L + z.',Nx,1);
  S.A1_xx = (2.0/L)*F;             S.A1_xz = -(1.0/L)*h.*FX;
  S.A2_xx = (1.0/L^2)*F.^2;        S.A2_xz = -(1.0/L^2)*h.*(F.*FX);
  S.A2_zz = (1.0/L^2)*(h.^2).*(FX.^2);
  S.B1_x = -(1.0/L)*FX;            S.B2_x = -(1.0/L^2)*F.*FX;
  S.B2_z = (1.0/L^2)*h.*(FX.^2);
  S.S1 = (2.0/L)*F;                S.S2 = (1.0/L^2)*F.^2;
  S.ell_bc = Nz+1; S.ell_if = 1; S.sgn_T = -1;
  S.phase = exp(-1i*gammap*z.');
end
% collocation matrices of solvebvp_colloc_fast(_lf), one per Fourier mode
% -- inverted once
gg = k2 - alphap.^2;
Ainv = zeros(Nz+1,Nz+1,Nx);
for j=1:Nx
  A = D2 + gg(j)*eye(Nz+1);
  if upper   % row 1: Robin (TBC) at z = L,  row end: Dirichlet at the interface z = 0
    A(end,:) = 0;  A(end,end) = 1;
    A(1,:) = D(1,:); A(1,1) = A(1,1) - 1i*gammap(j);
  else       % row 1: Dirichlet at the interface z = 0, row end: Robin (TBC) at z = -L
    A(end,:) = D(end,:); A(end,end) = A(end,end) + 1i*gammap(j);
    A(1,:) = 0;  A(1,1) = 1;
  end
  Ainv(:,:,j) = inv(A);
end
S.Ainv = Ainv;
S.T = T_dno(alpha,p,gamma,gammap,k2,Nx,M);   % p (= alphap for alpha = 0)
S.V = zeros(Nx,Nz+1,M+1,3);               % u_{n,m} for the last three n
S.tr_bc = zeros(Nx,M+1,N+1); S.tr_if = S.tr_bc; S.tz_if = S.tr_bc; S.D = S.tr_bc;
S.uhp = [];
end

function u_x = ddx(S,u)
u_x = ifft((1i*S.p).*fft(u));
end

function u_z = ddz(S,u)
u_z = u*S.DzLt;
end

function v = getV(S,m,n)
v = S.V(:,:,m+1,mod(n,3)+1);
end

function Fn = rhs(S,n,m)
% volume source of order (n,m): field_tfe_helmholtz_m_and_n(_lf).m, Fnm
Fn = zeros(S.Nx,S.Nz+1); al = S.alpha; g2 = S.g2;
if n>=1
  u1 = getV(S,m,n-1); ux = ddx(S,u1); uz = ddz(S,u1);
  Fn = Fn - ddx(S,S.A1_xx.*ux + S.A1_xz.*uz) - ddz(S,S.A1_xz.*ux) - S.B1_x.*ux ...
          - 2i*al*S.S1.*ux - g2*S.S1.*u1;
end
if m>=1
  u1 = getV(S,m-1,n);
  Fn = Fn - 2i*al*ddx(S,u1) - 2*g2*u1;
end
if n>=1 && m>=1
  u1 = getV(S,m-1,n-1);
  Fn = Fn - 2i*al*S.S1.*ddx(S,u1) - 2*g2*S.S1.*u1;
end
if n>=2
  u1 = getV(S,m,n-2); ux = ddx(S,u1); uz = ddz(S,u1);
  Fn = Fn - ddx(S,S.A2_xx.*ux + S.A2_xz.*uz) - ddz(S,S.A2_xz.*ux + S.A2_zz.*uz) ...
          - S.B2_x.*ux - S.B2_z.*uz - 2i*al*S.S2.*ux - g2*S.S2.*u1;
end
if m>=2
  Fn = Fn - g2*getV(S,m-2,n);
end
if n>=1 && m>=2
  Fn = Fn - g2*S.S1.*getV(S,m-2,n-1);
end
if n>=2 && m>=1
  u1 = getV(S,m-1,n-2);
  Fn = Fn - 2i*al*S.S2.*ddx(S,u1) - 2*g2*S.S2.*u1;
end
if n>=2 && m>=2
  Fn = Fn - g2*S.S2.*getV(S,m-2,n-2);
end
end

function G = dno(S,n,m,tz)
% dno_tfe_helmholtz_m_and_n(_lf).m at order (n,m), with tz = u_z at the interface
L = S.L; fx = S.fx_c; f = S.f;
if S.upper
  G = -tz;
  if n>=1
    G = G + fx.*ddx(S,S.tr_if(:,m+1,n)) + (1.0/L)*(f.*S.D(:,m+1,n));
  end
  if n>=2
    G = G - (1.0/L)*(f.*(fx.*ddx(S,S.tr_if(:,m+1,n-1)))) - fx.*(fx.*S.tz_if(:,m+1,n-1));
  end
else
  G = tz;
  if n>=1
    G = G - fx.*ddx(S,S.tr_if(:,m+1,n)) - (1.0/L)*(f.*S.D(:,m+1,n));
  end
  if n>=2
    G = G - (1.0/L)*(f.*(fx.*ddx(S,S.tr_if(:,m+1,n-1)))) + fx.*(fx.*S.tz_if(:,m+1,n-1));
  end
end
end

function [S,G] = partial(S,n,m)
% particular field of order (n,m) (zero interface data) and its DNO
if n==0 && m==0
  S.uhp = []; G = zeros(S.Nx,1); return;
end
Fn = rhs(S,n,m);
Jn = zeros(S.Nx,1);
for r=0:m-1
  Jn = Jn + S.sgn_T*ifft(S.T(:,m-r+1).*fft(S.tr_bc(:,r+1,n+1)));
end
if n>=1
  for r=0:m
    Jn = Jn - (1.0/S.L)*S.f.*ifft(S.T(:,m-r+1).*fft(S.tr_bc(:,r+1,n)));
  end
end
Fh = fft(Fn);                         % (Nx, Nz+1), Fourier in x
Fh(:,S.ell_bc) = fft(Jn);             % artificial-boundary (TBC) row
Fh(:,S.ell_if) = 0;                   % interface row: zero Dirichlet data
% uhp(j,:) = Ainv(:,:,j) * Fh(j,:).'  for all j at once (implicit expansion)
S.uhp = permute(sum(S.Ainv .* reshape(Fh.',1,S.Nz+1,S.Nx), 2), [3 1 2]);
tz = ifft(S.uhp*S.DzLt(:,S.ell_if));
G = dno(S,n,m,tz);
end

function S = complete(S,n,m,xi)
% add the homogeneous part exp(+-i gamma_p z) xi^ and store what later orders need
uh = S.phase .* fft(xi);
if ~isempty(S.uhp), uh = uh + S.uhp; end
u = ifft(uh);
S.V(:,:,m+1,mod(n,3)+1) = u;
S.tr_bc(:,m+1,n+1) = u(:,S.ell_bc);
S.tr_if(:,m+1,n+1) = u(:,S.ell_if);
S.tz_if(:,m+1,n+1) = u*S.DzLt(:,S.ell_if);
S.D(:,m+1,n+1) = dno(S,n,m,S.tz_if(:,m+1,n+1));
S.uhp = [];
end

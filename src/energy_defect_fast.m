function [ee,ru,rl] = energy_defect_fast(tau2,ubar_n_m,wbar_n_m,...
    d,alpha_bar,gamma_u_bar,gamma_w_bar,Eps,delta,...
    Nx,N,M,N_Eps,N_delta,SumType)
% energy_defect.m: Computes the energy defect (D) in order to inspect the conservation of energy.
%
%  Inputs:
%   tau2: a numerical constant representing TE or TM mode
%   ubar_n_m: a tensor representing approximate solution in upper layer
%   wbar_n_m: tensor representing approximate solution in lower layer
%   d: the periodicity of the grating
%   alpha_bar: a numerical constant
%   gamma_u_bar: a numerical constant in the upper field
%   gamma_w_bar: a numerical constant in the lower field
%   Eps: the physical error in the surface deformation
%   delta: the numerical error in the discretization of the frequency perturbation
%   N: the maximum number of Taylor orders for the interfacial perturbation
%   M: the maximum number of Taylor orders for the frequency perturbation
%   N_Eps: the maximum number of epsilon orders for the interfacial perturbation
%   N_Delta: the maximum number of delta orders for the frequency perturbation
%   SumType: boolean to control Taylor or Pade summation 
%
%  Outputs:
%   ee: an estimate of the total energy defect (which should approach 1)
%   ru: proportion of scattered energy in the upper field
%   rl: proportion of scattered energy in the lower field

if SumType ~= 1 && SumType ~= 2
  [ee,ru,rl] = energy_defect(tau2,ubar_n_m,wbar_n_m,d,alpha_bar,gamma_u_bar,gamma_w_bar,...
                             Eps,delta,Nx,N,M,N_Eps,N_delta,SumType);
  return;
end
Eps = Eps(:).'; delta = delta(:).';
[EE,DD] = ndgrid(Eps,delta);                  % (N_Eps, N_delta)
U = series_sum(ubar_n_m,EE(:).',DD(:).',N,M,SumType,Nx);   % (Nx, N_Eps*N_delta)
W = series_sum(wbar_n_m,EE(:).',DD(:).',N,M,SumType,Nx);
Uh = fft(U)/Nx; Wh = fft(W)/Nx;

ru = zeros(N_Eps,N_delta); rl = zeros(N_Eps,N_delta);
k_u_bar = sqrt(alpha_bar^2 + gamma_u_bar^2);
k_w_bar = sqrt(alpha_bar^2 + gamma_w_bar^2);
PropMode_u = zeros(Nx,1); PropMode_w = zeros(Nx,1);
for ell=1:N_delta
  s1 = 1 + delta(ell);
  [~,~,alpha_p,gamma_u_p] = setup_2d(Nx,d,s1*alpha_bar,s1*gamma_u_bar);
  [~,~,~,gamma_w_p] = setup_2d(Nx,d,s1*alpha_bar,s1*gamma_w_bar);
  gamma_i = gamma_u_p(1);
  % MATLAB compares complex numbers by their real parts
  PropMode_u(real(alpha_p.^2) < real((s1*k_u_bar)^2)) = 1;
  PropMode_w(real(alpha_p.^2) < real((s1*k_w_bar)^2)) = 1;
  cols = (ell-1)*N_Eps + (1:N_Eps);
  ru(:,ell) = sum((PropMode_u.*gamma_u_p./gamma_i) .* abs(Uh(:,cols)).^2, 1).';
  rl(:,ell) = tau2*sum((PropMode_w.*gamma_w_p./gamma_i) .* abs(Wh(:,cols)).^2, 1).';
end
ee = 1.0 - ru - rl;
end

% ------------------------------------------------------------------------------------------
function S = series_sum(c,E,Dl,N,M,SumType,Nx)
% c: (N+1, M+1, Nx) as passed to energy_defect (c(n+1,m+1,x) multiplies eps^n delta^m)
% returns S (Nx, P) = the series summed at the P points (E(k), Dl(k)) as padesum2/taylorsum_2_coeff
P = numel(E);
rho = sqrt(E.^2 + Dl.^2); theta = atan2(Dl,E);
K = min(N,M); Mp = floor(K/2);
cs = cos(theta); sn = sin(theta);
coeff = zeros(K+1,Nx,P);                      % ctilde_p for every x and point
for p=0:K
  Cp = zeros(Nx,p+1); Wp = zeros(p+1,P);
  for q=0:p
    Cp(:,q+1) = reshape(c(p-q+1,q+1,:),Nx,1);
    Wp(q+1,:) = cs.^(p-q).*sn.^q;
  end
  coeff(p+1,:,:) = reshape(Cp*Wp,1,Nx,P);
end
coeff = reshape(coeff,K+1,Nx*P);
r = reshape(repmat(rho,Nx,1),1,Nx*P);
if SumType == 1
  % taylorsum(coeff,rho,Mp): sum_{k=0}^{Mp} c_k rho^k
  S = horner(coeff(1:Mp+1,:),r);
elseif Mp == 0
  S = coeff(1,:);
else
  S = zeros(1,Nx*P);
  chunk = 20000;
  for i0 = 1:chunk:Nx*P
    idx = i0:min(Nx*P,i0+chunk-1);
    S(idx) = pade_batch(coeff(:,idx),r(idx),Mp);
  end
end
S = reshape(S,Nx,P);
end

function y = horner(a,x)
% sum_k a(k+1,:) .* x.^k  for every column
y = a(end,:);
for k=size(a,1)-1:-1:1
  y = y.*x + a(k,:);
end
end

function s = pade_batch(c,rho,Mp)
% padesum.m for many coefficient vectors c(:,j) at once
B = size(c,2);
H = zeros(Mp,Mp,B); rhs = zeros(Mp,B);
for i=1:Mp
  for j=1:Mp
    H(i,j,:) = reshape(c(Mp+i-j+1,:),1,1,B);     % toeplitz(c(M+1:2M), c(M+1:-1:2))
  end
  rhs(i,:) = -c(Mp+i+1,:);
end
bb = gauss_pp(H,rhs);                            % (Mp, B)
b = [ones(1,B); bb];
a = zeros(Mp+1,B);
for k=0:Mp
  for i=0:k
    a(k+1,:) = a(k+1,:) + b(i+1,:).*c(k-i+1,:);
  end
end
s = horner(a,rho)./horner(b,rho);
end

function x = gauss_pp(A,r)
% solve A(:,:,j) x(:,j) = r(:,j) for all j: Gaussian elimination with partial pivoting
n = size(A,1); B = size(A,3);
for k=1:n
  [~,piv] = max(abs(A(k:n,k,:)),[],1);
  piv = reshape(piv,1,B) + k - 1;
  for i=k+1:n
    sw = find(piv == i);
    if ~isempty(sw)
      tmp = A(k,:,sw); A(k,:,sw) = A(i,:,sw); A(i,:,sw) = tmp;
      tmp = r(k,sw); r(k,sw) = r(i,sw); r(i,sw) = tmp;
    end
  end
  for i=k+1:n
    l = A(i,k,:)./A(k,k,:);
    A(i,k:n,:) = A(i,k:n,:) - l.*A(k,k:n,:);
    r(i,:) = r(i,:) - reshape(l,1,B).*r(k,:);
  end
end
x = zeros(n,B);
for i=n:-1:1
  s = r(i,:);
  for j=i+1:n
    s = s - reshape(A(i,j,:),1,B).*x(j,:);
  end
  x(i,:) = s./reshape(A(i,i,:),1,B);
end
end

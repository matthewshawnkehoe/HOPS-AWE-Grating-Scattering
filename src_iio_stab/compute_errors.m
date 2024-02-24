function [relerr,nplot] = compute_errors(nu,Gn_tfe,Eps,N,Nx)

relerr = zeros(N+1,2);
nplot = zeros(N+1,1);

nu_tfe_taylor = Gn_tfe(:,0+1);
nu_tfe_pade = Gn_tfe(:,0+1);

nplot(0+1) = 0;
relerr(0+1,1) = norm(nu-nu_tfe_taylor,inf)/norm(nu,inf);
relerr(0+1,2) = norm(nu-nu_tfe_pade,inf)/norm(nu,inf);

for n=1:N
  M = floor(n/2);
  for j=1:Nx
    coeff = Gn_tfe(j,:).';
    nu_tfe_taylor(j) = taylorsum(coeff,Eps,n);
    nu_tfe_pade(j) = padesum(coeff,Eps,M);
  end
  nplot(n+1) = n;
  relerr(n+1,1) = norm(nu-nu_tfe_taylor,inf)/norm(nu,inf);
  relerr(n+1,2) = norm(nu-nu_tfe_pade,inf)/norm(nu,inf);
end

return;

%
% taylorsum
%

function [tsum] = taylorsum(c,Eps,N)
% taylorsum - Sums a truncated Taylor series.
%
% Inputs:
%
% c - Taylor series coefficients: [c_0,...,c_N]
% Eps - Value at which to sum
% N - Degree of truncated Taylor series (unused!)
%
% Outputs:
%
% tsum - Taylor sum evaluated at Eps
%
% DPN 2/8/18

tsum = polyval(c(N+1:-1:1),Eps);

return;

%
% padesum
%

function [psum,a,b] = padesum(c,Eps,M)
% padesum - Sums a truncated Taylor series via Pade approximation.
%
% Inputs:
%
% c - Taylor series coefficients: [c_0,...,c_N]
% Eps - Value at which to sum
% M - Numerator and denominator degree (M = N/2)
%
% Outputs:
%
% psum - Pade approximant evaluated at Eps
% a - Numerator coefficients: [a_0,...,a_M]
% b - Denominator coefficients: [1,b_1,...,b_M]
%
% Note: Assumes that N = 2*M (otherwise ignores c_N)
%
% DPN 2/8/18

if(M==0)
  a = c(1);
  b = [1];
else
  H = toeplitz(c(M+1:2*M-1+1),c(M+1:-1:1+1));
  bb = -H\c(M+1+1:2*M+1);
  b = [1;bb];
  aa = conv(b,c(0+1:M+1));
  a = aa(0+1:M+1);
end
psum = polyval(a(end:-1:1),Eps)./polyval(b(end:-1:1),Eps);

return;
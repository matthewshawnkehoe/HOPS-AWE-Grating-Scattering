function [U,W] = AInverse(Q,R,gammap,gammapw,Nx,tau2)
  
a = zeros(Nx,1);
b = zeros(Nx,1);

Q_hat = fft(Q);
R_hat = fft(R);

for j=1:Nx
  det_p = tau2*gammapw(j) + gammap(j);
  a(j) = ((tau2*gammapw(j))*Q_hat(j) + 1i*R_hat(j))/det_p;
  b(j) = ((-gammap(j))*Q_hat(j) + 1i*R_hat(j))/det_p;
end

U = ifft(a);
W = ifft(b);

return;
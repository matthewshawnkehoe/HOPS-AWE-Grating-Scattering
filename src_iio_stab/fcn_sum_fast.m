function [f] = fcn_sum_fast(SumType,f_n_m,Eps,delta,Nx,N,M)

f = zeros(Nx,1);

for j=1:Nx
  if(SumType==1)
    % Change to taylorsum_2_coeff if M,N >= 6
    f(j) = taylorsum_2_coeff(f_n_m(:,:,j),Eps,delta,N,M);
  elseif(SumType==2)
    f(j) = padesum2(f_n_m(:,:,j),Eps,delta,N,M);
  else 
    f(j) = padesum2_safe(f_n_m(:,:,j),Eps,delta,N,M);
  end
end

return;
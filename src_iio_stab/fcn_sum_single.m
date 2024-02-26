function [f] = fcn_sum_single(SumType,f_n,Eps,Nx,N)

f = zeros(Nx,1);

for j=1:Nx
  if(SumType==1)
    f(j) = taylorsum(f_n(:,j),Eps,N);
  elseif(SumType==2)
    f(j) = padesum(f_n(:,j),Eps,N);
  else 
    f(j) = padesum_safe(f_n(:,j),Eps,N);
  end
end

return;
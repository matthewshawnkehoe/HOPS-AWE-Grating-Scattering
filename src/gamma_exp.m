function [gamma_q_m] = gamma_exp(alpha_bar,alpha_bar_q,...
    gamma_bar,gamma_bar_q,k_bar,M)
% gamma_exp.m: Computes gamma_q_m through a series expansion of gamma_q.
%
%  Inputs:
%   alpha_bar: a numerical constant 
%   alpha_bar_q: a numerical constant generated for all wave numbers q
%   gamma_bar: a numerical constant 
%   gamma_bar_q: a numerical constant generated for all wave numbers q
%   k_bar: a numerical constant
%   M: the maximum number of Taylor orders for the frequency perturbation
%
%  Output:
%   gamma_q_m: the series expansion of gamma_q.

gamma_q_m = zeros(1,M+1);

gamma_q_m(0+1) = gamma_bar_q;

if(M>=1)
  gamma_q_m(1+1) = (k_bar^2-alpha_bar*alpha_bar_q)/gamma_q_m(0+1);
end

if(M>=2)
  gamma_q_m(2+1) = (gamma_bar^2-gamma_q_m(1+1)^2)/(2*gamma_q_m(0+1));
end

for m=3:M
  num_sum = 0;
  for r=1:m-1
    num_sum = num_sum - gamma_q_m(m-r+1)*gamma_q_m(r+1);
  end
  gamma_q_m(m+1) = num_sum/(2*gamma_q_m(0+1));
end

return;
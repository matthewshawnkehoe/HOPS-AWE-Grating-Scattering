function [] = plot_errors(fig_num,sum_type,N,M,...
    err_G,err_U,err_ubar,err_J,err_W,err_wbar)
% plot_errors.m: Plot errors for field solvers and DNOs G and J.
%
%  Inputs:
%   fig_num: figure number for plotting
%   sum_type: type of summation (Taylor or Padé)
%   N: the maximum number of Taylor orders for the interfacial perturbation
%   M: the maximum number of Taylor orders for the frequency perturbation
%   err_G: error for the upper layer DNO, G
%   err_U: error for the approximate solution (U) at the surface z=g(x) in 
%          the upper field
%   err_ubar: error for the approximate solution (u_bar) from the upper 
%             field solver at the top collocation point
%   err_J: error for the lower layer DNO, J
%   err_W: error for the approximate solution (W) at the surface z=g(x) in 
%          the lower field
%   err_wbar: error for the approximate solution (w_bar) from the lower 
%             field solver at the bottom collocation point
%
%  Outputs:
%   None

figure(fig_num);

subplot(2,3,1);
contourf([0:N],[0:M],log10(err_G).');
xlabel('$n$','interpreter','latex');
ylabel('$m$','interpreter','latex');
ss = sprintf('Error in $G$ (%s)',sum_type);
title(ss,'interpreter','latex');
colorbar;

subplot(2,3,2);
contourf([0:N],[0:M],log10(err_U).');
xlabel('$n$','interpreter','latex');
ylabel('$m$','interpreter','latex');
ss = sprintf('Error in $U$ (%s)',sum_type);
title(ss,'interpreter','latex');
colorbar;

subplot(2,3,3);
contourf([0:N],[0:M],log10(err_ubar).');
xlabel('$n$','interpreter','latex');
ylabel('$m$','interpreter','latex');
ss = sprintf('Error in $\\bar{u}$ (%s)',sum_type);
title(ss,'interpreter','latex');
colorbar;

subplot(2,3,4);
contourf([0:N],[0:M],log10(err_J).');
xlabel('$n$','interpreter','latex');
ylabel('$m$','interpreter','latex');
ss = sprintf('Error in $J$ (%s)',sum_type);
title(ss,'interpreter','latex');
colorbar;

subplot(2,3,5);
contourf([0:N],[0:M],log10(err_W).');
xlabel('$n$','interpreter','latex');
ylabel('$m$','interpreter','latex');
ss = sprintf('Error in $W$ (%s)',sum_type);
title(ss,'interpreter','latex');
colorbar;

subplot(2,3,6);
contourf([0:N],[0:M],log10(err_wbar).');
xlabel('$n$','interpreter','latex');
ylabel('$m$','interpreter','latex');
ss = sprintf('Error in $\\bar{w}$ (%s)',sum_type);
title(ss,'interpreter','latex');
colorbar;

return;
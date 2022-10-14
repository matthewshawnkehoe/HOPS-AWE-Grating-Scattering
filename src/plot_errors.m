function [] = plot_errors(fig_num,sum_type,N,M,...
    err_G,err_U,err_ubar,err_J,err_W,err_wbar)

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
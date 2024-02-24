% make_plots.m
%
% MATLAB script to plot relative error versus perturbation order.

clear all; close all;

SavePlots = 0;

% Interface: 1->U, 2->V_u, 3->V_ell, 4->W
Interface = 1;
% SumType: 1->Taylor, 2->Pade
SumType = 1;
M = SumType;

% Specify the four .mat files
files = {'three_10_0.005.mat', 'three_10_0.01.mat', 'three_10_0.05.mat', 'three_10_0.1.mat'};
num_files = length(files);

Eps_full = zeros(num_files, 1);
nplot = [];
relerr = [];

% Define line specifications
line_specs = {'b-o', 'g-*', 'r-<', 'c->', 'y-d', 'k-^'};

% Vectorized data loading and accumulation
for j = 1:num_files
    filename = files{j};
    load(filename);
    Eps_full(j) = Eps;

    if(Interface==1)
        nplot = [nplot nplot_U];
        relerr = [relerr relerr_U(:,M)];
    elseif(Interface==2)
        nplot = [nplot nplot_V_u];
        relerr = [relerr relerr_V_u(:,M)];
    elseif(Interface==3)
        nplot = [nplot nplot_V_ell];
        relerr = [relerr relerr_V_ell(:,M)];
    else
        nplot = [nplot nplot_W];
        relerr = [relerr relerr_W(:,M)];
    end
end

% First Plot: Relative Error versus N
figure(1);
for j = 1:num_files
    line_idx = mod(j-1, length(line_specs)) + 1;
    line_str = line_specs{line_idx};
    disp_name = sprintf('$\\varepsilon = %g$', Eps_full(j));
    
    semilogy(nplot(:, j), relerr(:, j), line_str, ...
        'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', disp_name);
    hold on;
end
title('Relative Error versus $N$','Interpreter','latex');
xlabel('$N$','Interpreter','latex');
ylabel('Relative Error','Interpreter','latex');
set(gca, 'fontsize', 16);
legend('show', 'Interpreter', 'latex', 'Location', 'northeast');

% Second Plot: Relative Error versus Epsilon
figure(2);
for n = 0:2:N
    line_idx = mod(n/2, length(line_specs)) + 1;
    line_str = line_specs{line_idx};
    disp_name = sprintf('$N = %d$', n);
    
    loglog(Eps_full, relerr(n+1, :), line_str, ...
        'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', disp_name);
    hold on;
end
title('Relative Error versus $\\\varepsilon$','interpreter','latex');
xlabel('$\\\varepsilon$','interpreter','latex');
ylabel('Relative Error','interpreter','latex');
set(gca, 'fontsize', 16);
legend('show', 'Interpreter', 'latex', 'Location', 'northeast');

if(SavePlots == 1)
    saveas(gcf, 'conv_N', 'epsc');
    saveas(gcf, 'conv_Eps', 'epsc');
end
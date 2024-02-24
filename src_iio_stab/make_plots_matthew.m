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

Eps_full = zeros(num_files,1);
nplot = [];
relerr = [];

% Define line specifications
line_specs = {'b-o', 'g-*', 'r-<', 'c->', 'y-d', 'k-^'};

for j=1:num_files
  filename = files{j};
  load(filename);
  Eps_full(j) = Eps;
  line_idx = mod(j-1, length(line_specs)) + 1;
  line_str = line_specs{line_idx};
  disp_name = sprintf('$\\varepsilon = %g$',Eps_full(j));

  if(Interface==1)
    nplot = [nplot nplot_U];
    relerr = [relerr relerr_U(:,M)];
    display(size(nplot));
    display(size(relerr));
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
  figure(1); hh1 = gca;
  semilogy(nplot(:,j),relerr(:,j),line_str,...
      'LineWidth',2,'MarkerSize',8,...
      'DisplayName',disp_name);
  axis([min(nplot(:,j)) max(nplot(:,j)) 1e-16 max(relerr(:,j))]);
  title('Relative Error versus $N$','Interpreter','latex');
  xlabel('$N$','Interpreter','latex');
  ylabel('Relative Error','Interpreter','latex');
  set(hh1,'fontsize',16);
  hold on;
end

for n=0:2:N
  line_idx = mod(n/2, length(line_specs)) + 1;
  line_str = line_specs{line_idx};
  disp_name = sprintf('$N = %d$',n);
  
  figure(2); hh2 = gca;
  loglog(Eps_full,relerr(n+1,:),line_str,...
      'LineWidth',2,'MarkerSize',8,...
      'DisplayName',disp_name);
  axis([min(Eps_full) max(Eps_full) 1e-16 max(relerr(:,j))])
  title('Relative Error versus $\\\varepsilon$','interpreter','latex');
  xlabel('$\\\varepsilon$','interpreter','latex');
  ylabel('Relative Error','interpreter','latex');
  set(hh2,'fontsize',16);
  hold on;
end

figure(1);
ll = legend('show','Location','northeast');
set(ll,'interpreter','latex');
figure(2);
ll = legend('show','Location','northeast');
set(ll,'interpreter','latex');

if(SavePlots==1)
  saveas(hh1,'conv_N','epsc');
  saveas(hh2,'conv_Eps','epsc');
end
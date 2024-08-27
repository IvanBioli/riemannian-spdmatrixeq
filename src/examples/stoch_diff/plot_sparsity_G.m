% Description: Stochastic Galerkin Matrix Equations
%   Plots the sparsity pattern of the matrices G

clear all; close all; clc;

% Using latex
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

TP = 5; load_stoch_diff;
%% Plot sparisty pattern G_i
fig = figure();
idx = 2:4;
t = tiledlayout(1, length(idx),'TileSpacing','Compact');
for i = idx
    nexttile;
    spy(G{i});
    txt = title("$G_{"+ num2str(i-1) + "}$"); txt.Interpreter= 'latex';
    xlabel('')
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
end

% Rescale figure
set(gcf, 'Units', 'centimeters');
textwidth = 14.68 * 1.3;
pos = get(gcf, 'Position'); pos(4) = 6; pos(3) = textwidth; set(gcf, 'Position', pos);
fontsize(gcf, 11, "points")
exportgraphics(fig, "../report/figures/sparsity_G.pdf")

%% Plot sparsity pattern Gmean and cholesky
fig = figure();
tiledlayout(1, 4,'TileSpacing','Compact');

nexttile;
spy(Gmean_old);
txt = title("$G$"); txt.Interpreter= 'latex';
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

nexttile;
spy(chol(Gmean_old));
txt = title("$\mathtt{chol}(G)$"); txt.Interpreter= 'latex';
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

nexttile;
spy(Gmean);
txt = title("$G$ \texttt{dissect} reord."); txt.Interpreter= 'latex';
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

nexttile;
spy(chol(Gmean));
txt = title("$\mathtt{chol}(G$ \texttt{dissect} reord.$)$"); txt.Interpreter= 'latex';
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

% Rescale figure
set(gcf, 'Units', 'centimeters');
textwidth = 14.68 * 1.3;
pos = get(gcf, 'Position'); pos(4) = 6; pos(3) = textwidth; set(gcf, 'Position', pos);
fontsize(gcf, 11, "points")
exportgraphics(fig, "../report/figures/sparsity_G_reordered.pdf")

%% Plot sparsity pattern reordered G_i
fig = figure();
idx = 2:4;
t = tiledlayout(1, length(idx),'TileSpacing','Compact');
for i = idx
    nexttile;
    spy(B{i});
    txt = title("$B_{"+ num2str(i-1) + "}$"); txt.Interpreter= 'latex';
    xlabel('')
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
end

% Rescale figure
set(gcf, 'Units', 'centimeters');
textwidth = 14.68 * 1.3;
pos = get(gcf, 'Position'); pos(4) = 6; pos(3) = textwidth; set(gcf, 'Position', pos);
fontsize(gcf, 11, "points")
exportgraphics(fig, "../report/figures/sparsity_B_reordered.pdf")
% Description: Produces plots for experiment "Finite difference discretization of 2D
% PDEs on square domain"

clearvars -except seed; clc;
%% Comparison of best Riemannian methods and CG with truncation
load results/manopt_PAPER_pdes_FD_n10000_a10_r12.mat;
names_manopt = names; results_manopt = results; manopt_fixed_rank = fixed_rank;
load results/cg_PAPER_pdesFD_n10000_a10.mat; results_cg = results; names_cg = names;
load results/cg_PAPER_cap12_pdesFD_n10000_a10.mat; results_cg_cap = results; names_cg_cap = names;
load results/manopt_PAPER_pdes_FD_n10000_a10_rank_adaptive.mat; results_adapt = results; names_adapt = names;
clearvars -except seed results_manopt names_manopt ...
    results_cg names_cg ...
    results_cg_cap names_cg_cap ...
    results_adapt names_adapt ...
    titletext manopt_fixed_rank savename
best_results = {};
best_names = {};
idx = 1;

% Get the best results from Riemannian optimization with fixed rank
best_manopt = zeros(1,38); best_manopt([20, 26, 28]) = 1;
for i = find(best_manopt)
    best_results{idx} = results_manopt{i};
    best_names{idx} = names_manopt{i};
    idx = idx + 1;
end

% Get the best results from CG with truncation
best_CG = [3];
for i = best_CG
    best_results{idx} = results_cg{i};
    best_names{idx} = names_cg{i};
    idx = idx + 1;
end
for i = best_CG
    best_results{idx} = results_cg_cap{i};
    best_names{idx} = names_cg_cap{i};
    idx = idx + 1;
end

% Get the best results from Riemannian optimization with adaptive rank
best_adapt = zeros(1,38); best_adapt(26) = 1;
for i = find(best_adapt)
    best_results{idx} = results_adapt{i};
    best_names{idx} = names_adapt{i};
    idx = idx + 1;
end

order = 1:length(best_names);
best_results = best_results(order);
best_names = best_names(order);

% Plot a comparison
plot_comparison_cg(best_results, best_names, manopt_fixed_rank, [], savename, true, true)
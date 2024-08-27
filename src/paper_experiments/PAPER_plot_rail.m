% Description: Produces plots for experiment "Modified Bilinear Rail problem"

clearvars -except k_rail; clc;
if ~exist("k_rail", "var")
    k_rail = 3;
end
%% Load results
load results/manopt_PAPER_rail_n5177_tangADIx8_r150.mat;
names_manopt = names; results_manopt = results; manopt_fixed_rank = fixed_rank;
load results/cg_PAPER_rail_n5177.mat; results_cg = results; names_cg = names;
load results/manopt_PAPER_rail_n5177_tangADIx8_rank_adaptive.mat; results_adapt = results; names_adapt = names;

clearvars -except k_rail ...
    results_manopt names_manopt ...
    names_cg results_cg ...
    names_adapt results_adapt ...
    titletext manopt_fixed_rank savename

%% Select results and plot
best_results = {};
best_names = {};
idx = 1;

% Get the best results from Riemannian optimization
best_manopt = [10, 13, 14];
for i = best_manopt
    best_results{idx} = results_manopt{i};
    best_names{idx} = names_manopt{i};
    idx = idx + 1;
end

% Get the best results from CG with truncation
best_CG = [1];
for i = best_CG
    best_results{idx} = results_cg{i};
    best_names{idx} = names_cg{i};
    idx = idx + 1;
end

% Get the best results from Riemannian optimization with adaptive rank
best_adapt = [14];
for i = best_adapt
    best_results{idx} = results_adapt{i};
    best_names{idx} = names_adapt{i};
    idx = idx + 1;
end

plot_comparison_cg(best_results, best_names, manopt_fixed_rank, [], savename, true, true)
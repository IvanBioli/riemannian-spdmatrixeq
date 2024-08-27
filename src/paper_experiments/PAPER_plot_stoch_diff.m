% Description: Produces plots for experiment "Stochastic Galerkin matrix equations"

clearvars -except seed; clc;
%% %%%%%%%%%%%%%%%%%%% TP5 %%%%%%%%%%%%%%%%%%%%%
% Comparison of best Riemannian methods and CG with truncation
load results/manopt_PAPER_stoch_diff_TP5_r55.mat;
names_manopt = names; results_manopt = results; manopt_fixed_rank = fixed_rank;
load results/cg_PAPER_stoch_diff_TP5.mat; results_cg = results; names_cg = names;
load results/MultiRB_TP5.mat; results_multirb = results; names_multirb = names;
load results/manopt_PAPER_stoch_diff_TP5_rank_adaptive.mat; results_adapt = results; names_adapt = names;
clearvars -except seed results_manopt names_manopt ...
    names_cg results_cg ...
    names_cg_cap results_cg_cap ...
    names_multirb results_multirb ...
    results_adapt names_adapt ...
    titletext manopt_fixed_rank savename
best_results = {};
best_names = {};
idx = 1;

% Get the best results from Riemannian optimization
best_manopt = [zeros(1,7), 0, 0, 0, 1, 0, 0, 0, zeros(1,10)];
for i = find(best_manopt)
    best_results{idx} = results_manopt{i};
    best_names{idx} = names_manopt{i};
    idx = idx + 1;
end

% Get the best results from CG with truncation
best_CG = [2];
for i = best_CG
    best_results{idx} = results_cg{i};
    best_names{idx} = names_cg{i};
    idx = idx + 1;
end

% Get the best results from MultiRB
best_MultiRB = [2];
for i = best_MultiRB
    best_results{idx} = results_multirb{i};
    best_names{idx} = names_multirb{i};
    idx = idx + 1;
end

% Get the best results from Riemannian optimization with adaptive rank
best_adapt = best_manopt;
for i = find(best_adapt)
    best_results{idx} = results_adapt{i};
    best_names{idx} = names_adapt{i};
    idx = idx + 1;
end

plot_comparison_cg(best_results, best_names, manopt_fixed_rank, [], savename, true, true)

%% %%%%%%%%%%%%%%%%%%% TP2 %%%%%%%%%%%%%%%%%%%%%
clearvars -except seed; clc;
% Comparison of best Riemannian methods and CG with truncation
load results/manopt_PAPER_stoch_diff_TP2_r180.mat;
names_manopt = names; results_manopt = results; manopt_fixed_rank = fixed_rank;
load results/cg_PAPER_stoch_diff_TP2.mat; results_cg = results; names_cg = names;
load results/MultiRB_TP2.mat; results_multirb = results; names_multirb = names;
load results/manopt_PAPER_stoch_diff_TP2_rank_adaptive.mat; results_adapt = results; names_adapt = names;
clearvars -except seed results_manopt names_manopt ...
    names_cg results_cg ...
    names_cg_cap results_cg_cap ...
    names_multirb results_multirb ...
    names_adapt results_adapt ...
    titletext manopt_fixed_rank savename
best_results = {};
best_names = {};
idx = 1;

% Get the best results from Riemannian optimization
best_manopt = [zeros(1,7), 0, 0, 0, 0, 1, 0, 0, zeros(1,10)];
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

% Get the best results from MultiRB
best_MultiRB = [2];
for i = best_MultiRB
    best_results{idx} = results_multirb{i};
    best_names{idx} = names_multirb{i};
    idx = idx + 1;
end

% Get the best results from Riemannian optimization with adaptive rank
best_adapt = [zeros(1,7), 0, 0, 0, 0, 1, 0, 0, zeros(1,10)];
for i = find(best_adapt)
    best_results{idx} = results_adapt{i};
    best_names{idx} = names_adapt{i};
    idx = idx + 1;
end

plot_comparison_cg(best_results, best_names, manopt_fixed_rank, [], savename, true, true)
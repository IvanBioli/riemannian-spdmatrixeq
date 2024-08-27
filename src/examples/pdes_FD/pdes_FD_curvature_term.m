% Description: Finite difference discretization of 2D PDEs with separable variables
%   Executes experiments about the eigenvalues of the curvature term and of
%   the Riemannian Hessian when optimizing with R-TR.
%   These experiments have not been included in the final version of the
%   thesis.

clear all; close all; clc;
%% Get all problem structs
debug = true;
whatruns_firstord = zeros(1,14);
whatruns_rtr = [1, 0, 0, 0, 0];
plot_and_save = false;  sanity_checks = false;
alpha = 10; n = 30; fixed_rank = 5;
manopt_pdes_FD;
% results = {results{29}, results{32}, results{34}, results{37}};
% names = {names{29}, names{32}, names{34}, names{37}};
results = {results{29}, results{34}};
names = {names{29} + " + tCG", names{34} + " + tCG"};

%% Riemannian Hessian with exact solver
opts = options_rtr;
opts.subproblemsolver =@trs_gep;
options_hess = struct();
options_hess.projehess = false;
problem.hess = @(X, eta, store) sylv_hess(X, eta, A, B, C, store, options_hess);
[~, ~, info, ~] = trustregions(problem, x0, opts);
results{5} = info;
names{5} = "Riemann. Hess + exact TRs";
plot_manopt_results(results, names, fixed_rank, [], [1, 1, 0, 0, 1], [], "dummy_2", false)

%% Projected Euclidean Hessian with exact solver
opts = options_rtr;
opts.subproblemsolver =@trs_gep;
options_hess = struct();
options_hess.projehess = true;
problem.hess = @(X, eta, store) sylv_hess(X, eta, A, B, C, store, options_hess);
[~, ~, info, ~] = trustregions(problem, x0, opts);
results{6} = info;
names{6} = "Proj. EHess + exact TRs";

%% Vectorized
% opts = options_rtr;
% opts.subproblemsolver =@trs_gep;
% problem.hess = @(X, eta, store) sylv_hess_vectorized(X, eta, A, B, C, store, problem.M);
% [~, ~, info, ~] = trustregions(problem, x0, opts);
% results{7} = info;
% names{7} = "Riemann. Hess w/o negative eigenvalues curv. + exact TRs";

%% Plots
save_flag = true;
savename = "n" + num2str(n) + "r" + num2str(fixed_rank);
if save_flag
    save("results/manopt_PDEs_curvature_" + savename + "_experiments.mat", ...
        "results", "names", "fixed_rank", "alpha", "n", "savename")
end

fig = figure();
for j = [1, 2, 5, 6]
    info = results{j};
    name = names{j};
    name = erase(name, "R-TR ");
    name = replace(name, "Riemann. Hess", "RHess");
    name = erase(name, ", $r=" + num2str(fixed_rank) + "$");
    lineopts =  {"DisplayName", name "Marker", "."}; %lineopts = lineopts_map(name, find(j == find(idx_j), 1), paperflag);
    semilogy([info.iter], [info.gradnorm], lineopts{:});
    hold on
end
xlabel('Iteration number $k$');
ylabel('Gradient norm');
set_figsize_legend(5, true, true)
if save_flag
    exportgraphics(fig, "../report/figures/manopt_PDEs_curvature_" + savename + "_experiments_1.pdf")
end

fig = figure();
for j = [5]
    info = results{j};
    idx = 5:length(info);
    name = names{j};
    lineopts =  {"Marker", "."}; %lineopts = lineopts_map(name, find(j == find(idx_j), 1), paperflag);
    plot([info(idx).iter], min([info(idx).eighess]), lineopts{:}, "DisplayName", "$\lambda_{\min}$ RHess", "Color", [0.9290    0.6940    0.1250]);
    hold on
    plot([info(idx).iter], min([info(idx).eigproj]), lineopts{:}, "DisplayName", "$\lambda_{\min}$ Proj. EHess", "Color", "black");
    plot([info(idx).iter], min([info(idx).eigcurv]), lineopts{:}, "DisplayName", "$\lambda_{\min}$ Curv. term", "Color", "magenta");
end
limsy=get(gca,'YLim');
set(gca,'Ylim',[limsy(1) limsy(2)+5]);
xlabel('Iteration number $k$');
set_figsize_legend(5, true, true)

% Save results
if save_flag
    exportgraphics(fig, "../report/figures/manopt_PDEs_curvature_" + savename + "_experiments_2.pdf")
end

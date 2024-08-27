% Description: Executes experiments for "Finite difference discretization of 2D PDEs
% on square domain"

clearvars -except seed;
alpha = 10; n = 10000;
%% CG with trunction
adi_iterations_vec = [8]; one_sided_prec = false; sanity_checks = false;
PAPER_cg_pdes_FD; close all; clearvars -except seed alpha n;
adi_iterations_vec = [8]; one_sided_prec = false; sanity_checks = false;
rank_cap = 12; sanity_checks = false; PAPER_cg_pdes_FD; close all; clearvars -except seed alpha n;

%% Riemannian optimization with fixed rank
whatruns_firstord_rgd = zeros(1,14);
whatruns_firstord_rcg = zeros(1,14); whatruns_firstord_rcg([6, 12, 14]) = 1;
whatruns_firstord = [whatruns_firstord_rgd, whatruns_firstord_rcg];
whatruns_rtr = zeros(1,10);
if exist("seed", "var") % Set seed if exists
    rng(seed);
end
sanity_checks = false; fixed_rank = 12; PAPER_manopt_pdes_FD; close all; clearvars -except seed alpha n ;

%% Riemannian optimization with adaptive rank
alpha = 10; n = 10000;

% Set rank adaptivity options
opts_adaptivity.r_start = 3;
opts_adaptivity.r_increase = 3;
tolres = 1e-6;
opts_adaptivity.tol_resnorm = tolres;
opts_adaptivity.maxiter_tot = 1000;
rel_tolgradnorm = opts_adaptivity.tol_resnorm / 2;
stopfun_opts.miniter = 5;
stopfun_opts.patience = 3;
stopfun_opts.gradnorm_decrease = 0;
stopfun_opts.resnorm_slopefactor = 0.75;
stopfun_opts.backward_window = 3;
opts_adaptivity.stopfun_options = stopfun_opts;

% Set solvers to test for rank adaptivity
whatruns_firstord = [zeros(1, 25), 1, 0, 0];
whatruns_rtr = zeros(1,10);
sanity_checks = false; rank_adaptivity = true;

% Run experiments with rank adaptivity
if exist("seed", "var") % Set seed if exists
    rng(seed);
end
PAPER_manopt_pdes_FD;

%% Plots
PAPER_plot_PDEs_FD;
clearvars -except seed ; close all;
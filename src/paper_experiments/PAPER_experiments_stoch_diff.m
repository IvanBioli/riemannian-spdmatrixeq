% Description: Executes experiments for "Stochastic Galerkin matrix equations"

%% TP 5
clearvars -except seed; close all; clc
TP = 5;

% CG with trunction
what_runs = [0, 1, 0, 0, 0];
PAPER_cg_stoch_diff; close all;
clearvars -except seed TP;

% MultiRB
tol = 5e-8;
p_option_vec = [2];
MultiRB_stoch_diff; close all;

%% Riemannian optimization with fixed rank
clearvars -except seed TP;
firstord_only = 2; rtr_only = 2;
whatruns_firstord = [0, 0, 0, 1, 0, 0, 0];
whatruns_rtr = zeros(1, 5);
fixed_rank = 55;
if exist("seed", "var") % Set seed if exists
    rng(seed);
end
PAPER_manopt_stoch_diff; close all;

%% Riemannian optimization with adaptive rank
clearvars -except seed TP whatruns_firstord whatruns_rtr firstord_only rtr_only;
% Set rank adaptivity options
rel_tolgradnorm = 1e-6 / 2;
opts_adaptivity.r_start = 5;
opts_adaptivity.r_increase = 10;
opts_adaptivity.tol_resnorm = 1e-6;
opts_adaptivity.maxiter_tot = 1000;
stopfun_opts.miniter = 3;
stopfun_opts.patience = 2;
stopfun_opts.gradnorm_decrease = 0;
stopfun_opts.resnorm_slopefactor = 0.75;
stopfun_opts.backward_window = 3;
opts_adaptivity.stopfun_options = stopfun_opts;
% Set solvers to test for rank adaptivity
firstord_only = 2; rtr_only = 2;
whatruns_firstord = [0, 0, 0, 1, 0, 0, 0];
whatruns_rtr = zeros(1, 5);
rank_adaptivity = true;
% Run experiments rank adaptivity
if exist("seed", "var") % Set seed if exists
    rng(seed);
end
PAPER_manopt_stoch_diff; %close all;

%% TP 2
clearvars -except seed; close all;
TP = 2;

% CG with trunction
what_runs = [0, 0, 1, 0, 0];
PAPER_cg_stoch_diff; close all;
clearvars -except seed TP;

% MultiRB
tol = 2e-7;
p_option_vec = [2];
MultiRB_stoch_diff;

% Riemannian optimization with fixed rank
clearvars -except seed TP;
whatruns_firstord = [0, 0, 0, 0, 1, 0, 0];
whatruns_rtr = zeros(1, 5);
firstord_only = 2; rtr_only = 2;
fixed_rank = 180;
if exist("seed", "var") % Set seed if exists
    rng(seed);
end
PAPER_manopt_stoch_diff; close all;

%% Riemannian optimization with adaptive rank
clearvars -except seed TP whatruns_firstord whatruns_rtr firstord_only rtr_only;
% Set rank adaptivity options
rel_tolgradnorm = 1e-5 / 2;
opts_adaptivity.r_start = 30;
opts_adaptivity.r_increase = 30;
opts_adaptivity.tol_resnorm = 2*1e-5;
opts_adaptivity.maxiter_tot = 1000;
stopfun_opts.miniter = 3;
stopfun_opts.patience = 2;
stopfun_opts.gradnorm_decrease = 0;
stopfun_opts.resnorm_slopefactor = 0.75;
stopfun_opts.backward_window = 3;
opts_adaptivity.stopfun_options = stopfun_opts;
% Run experiments rank adaptivity
rank_adaptivity = true;
if exist("seed", "var") % Set seed if exists
    rng(seed);
end
PAPER_manopt_stoch_diff;

%% Plots
clearvars -except seed;
PAPER_plot_stoch_diff;

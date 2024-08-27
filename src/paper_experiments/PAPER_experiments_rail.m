% Description: Executes experiments for "Modified Bilinear Rail problem"

%% CG with trunction
clear all -except seed; close all; clc;
k_rail = 3; adi_iterations_vec_two = [16]; PAPER_cg_rail;

%% Riemannian optimization with fixed rank
clearvars -except k_rail seed;
fixed_rank = 150;

k_tangADI = 8;
whatruns_firstord = zeros(1, 14); whatruns_firstord([10, 13, 14]) = 1;
whatruns_rtr = zeros(1,10);
PAPER_manopt_rail;

%% Riemannian optimization with adaptive rank
clearvars -except k_rail seed;

% Set rank adaptivity options
opts_adaptivity.r_start = 50;
opts_adaptivity.r_increase = 25;
opts_adaptivity.tol_resnorm = 1e-6;
opts_adaptivity.maxiter_tot = 1000;
tolres = opts_adaptivity.tol_resnorm;
rel_tolgradnorm = opts_adaptivity.tol_resnorm / 1.5;
stopfun_opts.miniter = 4;
stopfun_opts.patience = 2;
stopfun_opts.gradnorm_decrease = 0;
stopfun_opts.resnorm_slopefactor = 0.75;
stopfun_opts.backward_window = 3;

opts_adaptivity.stopfun_options = stopfun_opts;

% Set solvers to test for rank adaptivity
k_tangADI = 8;
whatruns_firstord = [zeros(1,13), 1];
whatruns_rtr = zeros(1,10);

% Run experiments with rank adaptivity
rank_adaptivity = true;
PAPER_manopt_rail;
PAPER_plot_rail;

%% Plots
PAPER_plot_rail;
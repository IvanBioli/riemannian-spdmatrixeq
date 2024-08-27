% Description: Stochastic Galerkin Matrix Equations
%   Executes experiments with CG with truncation

clearvars -except seed TP what_runs rank_cap; clc;
%% GENERATING PROBLEM INSTANCE
% Load data
% TP = 2;
load_stoch_diff;

%% DEFINE CG OPTIONS
% Define the operator
A_op =@(X) sylv_op_lr(A, B, X);

% Define PCG options
options.maxiter = 35;
if TP == 2
    rel_tol = 1e-5;
else
    rel_tol = 1e-6;
end
options.tol = rel_tol;
options.truncate_R = true;
options.truncate_Q = false;
options.maxrank = 500;

% Define the truncation operator
opts_trunc_X.method = ["relF"]; opts_trunc_X.tol = [0.0025 * rel_tol];
opts_trunc_R.method = ["relF", "absF_fact"]; opts_trunc_R.tol = [0.1 * rel_tol, 0.001 * rel_tol];
opts_trunc_Q.method = ["relF", "absF_fact"]; opts_trunc_Q.tol = [0.1 * rel_tol, 0.001 * rel_tol];
opts_trunc_P.method = ["relF", "absF_fact"]; opts_trunc_P.tol = [0.1 * rel_tol, 0.001 * rel_tol];

if exist("rank_cap", "var")
    opts_trunc_X.method = [opts_trunc_X.method, "rank"];
    opts_trunc_X.tol = [opts_trunc_X.tol, rank_cap];
    opts_trunc_R.method = [opts_trunc_R.method, "rank"];
    opts_trunc_R.tol = [opts_trunc_R.tol, rank_cap];
    opts_trunc_Q.method = [opts_trunc_Q.method, "rank"];
    opts_trunc_Q.tol = [opts_trunc_Q.tol, rank_cap];
    opts_trunc_P.method = [opts_trunc_P.method, "rank"];
    opts_trunc_P.tol = [opts_trunc_P.tol, rank_cap];
    rank_cap_str = "cap" + num2str(rank_cap) + "_";
    cap_legend = " cap $r = "+ num2str(rank_cap) + "$ +";
else
    rank_cap_str = "";
    cap_legend = "";
end
options.opts_trunc_X = opts_trunc_X;
options.opts_trunc_R = opts_trunc_R;
options.opts_trunc_Q = opts_trunc_Q;
options.opts_trunc_P = opts_trunc_P;
T_op = @(X, opts) truncation(X, opts);

%% RUN (P)CG
if ~exist("what_runs", "var")
    what_runs = [0, 1, 1, 0, 0];
end
% Run CG
if what_runs(1)
    M_inv_op =@(X) X;
    [~, info] = pcg_trunc_lr(A_op, M_inv_op, C, T_op, options);
    results{1} = info;
    names{1} = "CG";
end

% Run pCG with truncation + preconditioner P(X) = K0*X
if what_runs(2)
    M_inv_op =@(X) meanbased_precond(X, Amean_d);
    [~, info] = pcg_trunc_lr(A_op, M_inv_op, C, T_op, options);
    info = add_time(info, time_precomp_A);
    results{2} = info;
    names{2} = "CG trunc. +" + cap_legend + " $\mathcal{P}^{(1)}X = K_0X$";
end

% Run pCG with truncation + preconditioner P(X) = K*X*G
if what_runs(3)
    M_inv_op =@(X) meanbased_precond_stoch(X, Amean_d, Gmean_d);
    [~, info] = pcg_trunc_lr(A_op, M_inv_op, C, T_op, options);
    info = add_time(info, time_precomp_A + time_precomp_G);
    results{3} = info;
    names{3} = "CG trunc. +" + cap_legend + " $\mathcal{P}^{(2)}X = K_0XG$";
end

% Run pCG (full rank) + preconditioner P(X) = K_0*X
if what_runs(4)
    options_fullrank = struct('maxiter', options.maxiter, 'tol', options.tol, ...
        'lowrank', false);
    T_op_fullrank =@frobnorm;
    A_op_fullrank =@(X) sylv_op(A, B, X);
    M_inv_op_fullrank =@(X) meanbased_precond(X, Amean_d);
    [~, info] = pcg_trunc_lr(A_op_fullrank, M_inv_op_fullrank, ...
        full(C.L * C.R'), T_op_fullrank, options_fullrank);
    info = add_time(info, time_precomp_A);
    results{4} = info;
    names{4} = "CG (full rank) + $\mathcal{P}^{(1)}X = K_0X$";
end

% Run pCG (full rank) + preconditioner P(X) = K*X*G
if what_runs(5)
    options_fullrank = struct('maxiter', options.maxiter, 'tol', options.tol, ...
        'lowrank', false);
    T_op_fullrank =@frobnorm;
    A_op_fullrank =@(X) sylv_op(A, B, X);
    M_inv_op_fullrank =@(X) meanbased_precond_stoch(X, Amean_d, Gmean_d);
    [~, info] = pcg_trunc_lr(A_op_fullrank, M_inv_op_fullrank, ...
        full(C.L * C.R'), T_op_fullrank, options_fullrank);
    info = add_time(info, time_precomp_A + time_precomp_G);
    results{5} = info;
    names{5} = "CG (full rank) + $\mathcal{P}^{(2)}X = K_0XG$";
end

%% PLOTS
titletext = "Test Problem " + num2str(TP);
savename = rank_cap_str + "stoch_diff_TP" + num2str(TP);
plot_cg_results(results, names, titletext, savename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = meanbased_precond(X, Amean_d)
% meanbased_precond - Applies to X the preconditioner P(X) = K * X
%
% Syntax:
%   [X] = meanbased_precond(X, Amean_d)
%
% Inputs:
%   - X      : Matrix to which the preconditioner is applied. Can be a full
%              matrix or a structure with fields U and V representing the
%              matrix X.U * X.V'
%   - Amean_d: Decomposition of the matrix K


if isstruct(X)
    X.U = Amean_d \ X.U;
else
    X = Amean_d \ X;
end
end

function X = meanbased_precond_stoch(X, Amean_d, Gmean_d)
% meanbased_precond_stoch - Applies to X the preconditioner P(X) = K*X*G
%
% Syntax:
%   [X] = meanbased_precond_stoch(X, Amean_d, Gmean_d)
%
% Inputs:
%   - X      : Matrix to which the preconditioner is applied. Can be a full
%              matrix or a structure with fields U and V representing the
%              matrix X.U * X.V'
%   - Amean_d: Decomposition of the matrix K
%   - Gmean_d: Decomposition of the matrix G


if isstruct(X)
    X.U = Amean_d \ X.U;
    X.V = Gmean_d \ X.V;
else
    X = Amean_d \ X / Gmean_d;
end
end

function [X, r] = frobnorm(X, ~)
% frobnorm - Computes the frobenius norm of a full matrix X, and returns X,
% with a sintax compatible with the truncation operator of pcg_trunc_lr
%
% Syntax:
%   [X, r] = frobnorm(X, ~)

r = norm(X, 'fro');
end

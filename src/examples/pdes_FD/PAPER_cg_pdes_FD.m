% Description: Finite difference discretization of 2D PDEs with separable variables
%   Executes experiments with CG with truncation

clearvars -except seed alpha n sanity_checks adi_iterations_vec ...
    rank_cap one_sided_prec; clc;
%% GENERATING PROBLEM INSTANCE
% Parameters of the problem
% alpha = 10;
n_terms = 3;
% n = 1000;

% Assembling the system
[data, u] = example1_paper(alpha, n_terms);
[A, B, C] = assemble_pdes_FD(n+1, data);

% Whether to run the full-rank version or not
% sanity_checks = false;
% one_sided_prec = false;
%% DEFINE CG OPTIONS
% Define the operator
A_op =@(X) sylv_op_lr(A, B, X);

% Precomputation for ADI, P(X) = TX+XT preconditioner
e = ones(n, 1); T = data.k_mean * spdiags([-e 2*e -e], -1:1, n, n);
a = data.k_mean * 2 * (1 - cos(pi / (n+1))); % = eigs(T, 1, 'smallestreal');
b = data.k_mean * 2 * (1 - cos(n * pi / (n+1))); % = eigs(T, 1, 'largestreal');

% Precomputations for ADI, P(X) = AXD+EXB preconditioner
[A_dom, B_dom, ~] = assemble_pdes_FD(n+1, data.dominant);
A_ = A_dom{3}; D_ = B_dom{3}; E_ = A_dom{4}; B_ = B_dom{4};
E_d = decomposition(E_, 'diagonal'); D_d = decomposition(D_, 'diagonal');
time_shifts = tic();
b_dom = eigs(A_, E_, 1, 'largestabs', 'Tolerance', 1e-4);
a_dom = eigs(A_, E_, 1, 'smallestabs', 'Tolerance', 1e-4);
d_dom = eigs(B_, D_, 1, 'largestabs', 'Tolerance', 1e-4);
c_dom = eigs(B_, D_, 1, 'smallestabs', 'Tolerance', 1e-4);
time_shifts = toc(time_shifts);

% Define PCG options
rel_tol = 1e-6;
options.maxiter = 100;
options.tol = rel_tol;
options.truncate_R = true;
options.truncate_Q = false;
options.truncate_P = true;

% Define the truncation operators
opts_trunc_X.method = ["relF"]; opts_trunc_X.tol = [0.01 * rel_tol];
opts_trunc_R.method = ["relF", "absF_fact"]; opts_trunc_R.tol = [0.1 * rel_tol, 0.001 * rel_tol];
opts_trunc_Q.method = ["relF", "absF_fact"]; opts_trunc_Q.tol = [0.1 * rel_tol, 0.001 * rel_tol];
opts_trunc_P.method = ["relF", "absF_fact"]; opts_trunc_P.tol = [0.1 * rel_tol, 0.001 * rel_tol];

if exist("rank_cap", "var")
    opts_trunc_X.method = [opts_trunc_X.method, "rank"];
    opts_trunc_X.tol = [opts_trunc_X.tol, rank_cap];
    opts_trunc_R.method = [opts_trunc_R.method, "rank"];
    opts_trunc_R.tol = [opts_trunc_R.tol];
    opts_trunc_Q.method = [opts_trunc_Q.method, "rank"];
    opts_trunc_Q.tol = [opts_trunc_Q.tol];
    opts_trunc_P.method = [opts_trunc_P.method, "rank"];
    opts_trunc_P.tol = [opts_trunc_P.tol];
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

%% RUN CG
if ~exist("adi_iterations_vec", "var")
    adi_iterations_vec = [0, 1, 2, 4, 8, 16];
end
if adi_iterations_vec(1) == 0
    % No preconditioner
    M_inv_op = @(X) X;
    [~, info] = pcg_trunc_lr(A_op, M_inv_op, C, T_op, options);
    results{1} = info;
    names{1} = "CG with trunc. (no preconditioner)";
end

% ADI preconditioner
for i=1:length(adi_iterations_vec)
    adi_niter = adi_iterations_vec(i);
    
    if one_sided_prec
        % Run pCG with preconditioner P(X) = TX+XT
        [p, q] = zolotarev_poles(adi_niter, a, b, a, b);
        M_inv_op = @(X) fadi_precond(X, T, T, p, q);
        [~, info] = pcg_trunc_lr(A_op, M_inv_op, C, T_op, options);
        results{2*i} = info;
        names{2*i} = "CG trunc. +" + cap_legend + " $\mathcal{P}^{(1)} X =TX+XT$ fADI$\times " + num2str(adi_niter) + "$";
    end
    
    % Run pCG with preconditioner P(X) = AXD+EXB
    [p, q] = zolotarev_poles(adi_niter, a_dom, b_dom, c_dom, d_dom);
    M_inv_op = @(X) fadi_precond(X, A_, B_, p, q, D_, E_);
    [~, info] = pcg_trunc_lr(A_op, M_inv_op, C, T_op, options);
    info = add_time(info, time_shifts);
    results{2*i+1} = info;
    names{2*i+1} = "CG trunc. +" + cap_legend + " $\mathcal{P}^{(2)} X =AXD+EXB$ fADI$\times " + num2str(adi_niter) + "$";
end

if sanity_checks
    % Run full rank pCG (no trucation) with P(X) = TX + XT and exact Sylvester
    options_fullrank = struct('maxiter', options.maxiter, 'tol', options.tol, ...
        'lowrank', false);
    T_op_fullrank =@(X) X;
    A_op_fullrank =@(X) sylv_op(A, B, X);
    M_inv_op_fullrank =@(X) sylvester(full(T), full(T), X);
    [~, info_fullrank] = pcg_trunc_lr(A_op_fullrank, M_inv_op_fullrank, ...
        full(C.L * C.R'), T_op_fullrank, options_fullrank);
    results{length(results)+1} = info_fullrank;
    names{length(names)+1} = "CG (full rank) + $\mathcal{P}^{(1)} X = TX+XT$ exact";
    
    % Run full rank pCG (no trucation) with P(X) = AXB + EXD and exact
    % Sylvester
    options_fullrank = struct('maxiter', options.maxiter, 'tol', options.tol, ...
        'lowrank', false);
    T_op_fullrank =@(X) X;
    A_op_fullrank =@(X) sylv_op(A, B, X);
    M_inv_op_fullrank =@(X) bartelsStewart(A_, D_, E_, B_, X, false, false);
    [~, info_fullrank] = pcg_trunc_lr(A_op_fullrank, M_inv_op_fullrank, ...
        full(C.L * C.R'), T_op_fullrank, options_fullrank);
    results{length(results)+1} = info_fullrank;
    names{length(names)+1} = "CG (full rank) + $\mathcal{P}^{(2)} X = AXD+EXB$ (exact)";
end

%% PLOTS
titletext = "$n = " + n + "$, $a =" + alpha + "$";
savename = "PAPER_" + rank_cap_str + "pdesFD_n" + num2str(n) + "_a" + num2str(alpha);
plot_cg_results(results, names, titletext, savename)
% Description: Reachability Gramian of a bilinear control system
%   Executes experiments with CG with truncation

clearvars -except seed k_rail no_precond adi_iterations_vec_one ...
    adi_iterations_vec_two; clc;
%% GENERATING PROBLEM INSTANCE
% k_rail = 4;
load_rail;

%% DEFINE CG OPTIONS
% Define the operator
A_op =@(X) lyap_op_lr(A, X, M, N);

% Precomputations for the shifts
I = speye(n);
time_shifts_onesided = tic();
a_onesided = eigs(A, 1, 'smallestabs');
b_onesided = eigs(A, 1, 'largestabs');
time_shifts_onesided = toc(time_shifts_onesided);

time_shifts_twosided = tic();
a_twosided = eigs(A, M, 1, 'smallestabs');
b_twosided = eigs(A, M, 1, 'largestabs');
time_shifts_twosided = toc(time_shifts_twosided);

% Define PCG options
rel_tol = 1e-6;
options.maxrank = 1000;
options.maxiter = 50;
options.tol = rel_tol;
options.truncate_R = true;
options.truncate_Q = false;
options.truncate_P = true;

% Define the truncation operators
opts_trunc_X.method = ["relF"]; opts_trunc_X.tol = [0.00025 * rel_tol];
opts_trunc_R.method = ["relF", "absF_fact"]; opts_trunc_R.tol = [0.01 * rel_tol, 0.0001 * rel_tol];
opts_trunc_P.method = ["relF", "absF_fact"]; opts_trunc_P.tol = [0.01 * rel_tol, 0.0001 * rel_tol];
options.opts_trunc_X = opts_trunc_X;
options.opts_trunc_R = opts_trunc_R;
options.opts_trunc_P = opts_trunc_P;
T_op = @(X, opts) truncation(X, opts);

%% RUN CG
if ~exist("no_precond", "var")
    no_precond = false;
end
if ~exist("adi_iterations_vec_one", "var")
    adi_iterations_vec_one = [];
end
if ~exist("adi_iterations_vec_two", "var")
    adi_iterations_vec_two = [16];
end

% No preconditioner
if no_precond
    M_inv_op = @(X) X;
    [~, info] = pcg_trunc_lr(A_op, M_inv_op, B, T_op, options);
    results{1} = info;
    names{1} = "CG with trunc. (no preconditioner)";
end
idx_shift = no_precond;

% Run pCG with preconditioner P(X) = AX+XA and fADI
for i=1:length(adi_iterations_vec_one)
    adi_niter = adi_iterations_vec_one(i);
    p = zolotarev_poles(adi_niter, a_onesided, b_onesided);
    M_inv_op =@(X) adi_lyap(A, X, p);
    
    [~, info] = pcg_trunc_lr(A_op, M_inv_op, B, T_op, options);
    info = add_time(info, time_shifts_onesided + time_precomp);
    
    results{idx_shift+i} = info;
    names{idx_shift+i} = "CG trunc. + $\mathcal{P}^{(1)} X =AX+XA$ fADI$\times " + num2str(adi_niter) + "$";
end
idx_shift = idx_shift + length(adi_iterations_vec_one);

% Run pCG with preconditioner P(X) = AXM+MXA and fADI
for i=1:length(adi_iterations_vec_two)
    adi_niter = adi_iterations_vec_two(i);
    p = zolotarev_poles(adi_niter, a_twosided, b_twosided);
    M_inv_op =@(X) adi_lyap(A, X, p, struct(), M);
    
    [X, info] = pcg_trunc_lr(A_op, M_inv_op, B, T_op, options);
    info = add_time(info, time_shifts_twosided + time_precomp);
    
    results{idx_shift+i} = info;
    names{idx_shift+i} = "CG trunc. + $\mathcal{P}^{(2)} X =AXM+MXA$ fADI$\times " + num2str(adi_niter) + "$";
end

%% PLOTS
titletext = "Rail profile problem, $n = " + n + "$";
savename = "PAPER_rail_n" + num2str(n);
plot_cg_results(results, names, titletext, savename)

%%
% X_ = X.V * X.D * X.V';
% Lyap_res = A * X_ * M + M * X_ * A - B*B';
% norm(Lyap_res, 'fro') / norm(B*B', 'fro')
function [stats, store] = sylv_statsfun(problem, x, stats, store)
% sylv_statsfun - Routine in charge of logging additional
% statistics in Manopt's info struct-array, when solving the generalized
% Sylvester equation by optimization of the manifold of fixed-rank
% matrices. Time spent here is not counted in the solver's time. The
% function's signature is compatible with Manopt.
%   Logs: residual norm (res_norm), relative residual norm(res_norm_rel),
%   minimum/maximum singular value of the current iterate
%   (sigma_min/sigma_max), and if debug also the eigenvalues of the
%   Riemannian Hessian, of the curvature term and of the projected
%   Euclidean Hessian

if isfield(problem, 'debug') && problem.debug
    debug = true;
else
    debug = false;
end
stats.sigma_min = min(diag(x.S));
stats.sigma_max = max(diag(x.S));
if isstruct(x)
    [m, r] = size(x.U); n = size(x.V, 1); ell = size(store.AU, 3);
    residual.U = [reshape(bsxfun(@times, store.AU, diag(x.S)'), [m, r*ell]), -problem.C.L];
    residual.V = [reshape(store.BV, [n, r*ell]), problem.C.R];
    stats.res_norm = stable_norm_fact(residual);
    start_normest = tic();
    stats.res_norm_est = hutch_normest(residual, 5);
    time_normest = toc(start_normest);
    stats.time = stats.time + time_normest;
else % Full rank case, used only for sanity checks
    stats.res_norm = norm(store.grad__, 'fro');
    stats.res_norm_est = stats.res_norm_est(store.grad__, 5);
end
stats.res_norm_slope = nan; % Dummy value to be replaced
if isfield(problem, 'rhs_norm')
    stats.res_norm_rel = stats.res_norm / problem.rhs_norm;
    stats.res_norm_rel_est = stats.res_norm_est / problem.rhs_norm;
end
if debug
    [~, store] = sylv_hess_vectorized(x, problem.M.randvec(x), problem.A, problem.B, problem.C, store, problem.M);
    stats.eigproj = store.eigproj;
    stats.eigcurv = store.eigcurv;
    stats.eighess = store.eighess;
end
end


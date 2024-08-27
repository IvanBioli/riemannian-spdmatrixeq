function stats = lyap_psd_statsfun(problem, X, stats, store)
% lyap_psd_statsfun - Routine in charge of logging additional
% statistics in Manopt's info struct-array, when solving the generalized
% Lyapunov equation by optimization of the manifold of fixed-rank positive
% semidefinite matrices. Time spent here is not counted in the solver's
% time. The function's signature is compatible with Manopt.
%   Logs: residual norm (res_norm), relative residual norm(res_norm_rel)

    % Assemble L(X) - C in factored form
    [n, r] = size(X.V); ell = length(problem.N); B = problem.B;
    residual = struct();
    residual.V = [(store.AV + store.MV) / sqrt(2), ...
        (store.AV - store.MV) / sqrt(2), ...
        reshape(store.NV, [n, r*ell]), ...
        B];
    residual.D = blkdiag(kron(...
        spdiags([1;-1;-ones(ell, 1)], 0, ell+2, ell+2), ...
        X.D), ...
        -speye(size(B,2)));
    stats.res_norm = stable_norm_fact(residual);
    start_normest = tic();
    stats.res_norm_est = hutch_normest(residual, 5);
    time_normest = toc(start_normest);
    stats.time = stats.time + time_normest;

    stats.res_norm_slope = nan; % Dummy value to be replaced
    if isfield(problem, 'rhs_norm')
        stats.res_norm_rel = stats.res_norm / problem.rhs_norm;
        stats.res_norm_rel_est = stats.res_norm_est / problem.rhs_norm;
    end
end


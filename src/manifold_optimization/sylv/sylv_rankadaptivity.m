function [X, info_tot] = sylv_rankadaptivity(problem, solver, M_handle, options, D_d, E_d)
% sylv_rankadaptivity - Riemannian Rank Adaptive Algorithm
% obtained alternating between fixed-rank Riemannian optimization and
% rank-update steps. The fixed-rank manifold is the matrix of fixed-rank
% matrices, possibly with the tools in the metric induced by P(X) = EXD.
%
% Syntax:
%   [X, info_tot] = sylv_rankadaptivity(problem, solver, M_handle, options, D_d, E_d)
%
% Inputs:
%   - problem : Problem structure
%   - solver  : Solver for fixed-rank Riemannian optimization
%   - M_handle: Handle for defining the manifold given the rank
%   - options : Options structure. See mergeDefaultOptions.
%   - D_d     : Decomposition of the matrix D
%   - E_d     : Decomposition of the matrix E
%
% Outputs:
%   - X       : Final solution
%   - info_tot: Info structure
%

% Get the options
options = mergeDefaultOptions(options);
r_start = options.r_start;
r_increase = options.r_increase;
opts_solver = options.opts_solver;
tol_resnorm = options.tol_resnorm;
maxiter_tot = options.maxiter_tot;
stopfun_options = options.stopfun_options;
epsilon_sigma = stopfun_options.epsilon_sigma;

% Create global variable for heuristic stopping criterion
global slope_history
slope_history = zeros(1, 1000);

% Initialize the manifold with rank r_start
r_curr = r_start; err_frob = 1;
if ~isfield(options, "X0")
    problem.M = M_handle(r_start);
    X0 = problem.M.rand(); % Random starting point
else
    X0 = options.X0;
end
info_tot = []; % Info array struct
metric_EXD = isfield(problem.M, "qr_E"); % Get if we are using the metric induced by EXD
% Main loop
iter = 0;
while err_frob > tol_resnorm && iter < maxiter_tot
    % Define the manifold with rank r_curr
    problem.M = M_handle(r_curr);
    problem.linesearch = @(X, d, store) sylv_linesearch_init(X, d, store, problem);
    
    % Solve the problem on matrices of rank r_curr
    update_opts_solver(); flag_stopfun = false;
    [X, ~, info, ~, store] = solver(problem, X0, opts_solver);
    aux = num2cell(r_curr * ones(size(info)));
    [info.rank] = aux{:};
    
    % Update the rank
    s = diag(X.S);
    if flag_stopfun && trunc_rank(s) < r_curr % If stopped because iterate is rank deficient
        r_curr = trunc_rank(s);
        X0.U = X.U(:, 1:r_curr);
        X0.V = X.V(:, 1:r_curr);
        X0.S = X.S(1:r_curr, 1:r_curr);
        if metric_EXD
            X0.EU = X.EU(:, 1:r_curr); X0.DV = X.DV(:, 1:r_curr);
        end
        time_update = 0;
    else
        % Assemble the factored residual from already computed quantities
        [m, r] = size(X.U); n = size(X.V, 1); ell = size(store.AU, 3);
        residual.U = [reshape(bsxfun(@times, store.AU, sqrt(diag(X.S))'), [m, r*ell]), -problem.C.L];
        residual.V = [reshape(bsxfun(@times, store.BV, sqrt(diag(X.S))'), [n, r*ell]), problem.C.R];
        timetic = tic;
        
        % Euclidean gradient
        egrad = residual;
        if metric_EXD
            egrad.U = E_d \ egrad.U;
            egrad.V = D_d \ egrad.V;
        end
        
        % Component of the euclidean gradient orthogonal to the tangent
        perp_egrad.U = egrad.U - X.U * [reshape(bsxfun(@times, store.UtAU, sqrt(diag(X.S))'), [r, r*ell]), -store.CLtU'];
        perp_egrad.V = egrad.V - X.V * [reshape(bsxfun(@times, store.VtBV, sqrt(diag(X.S))'), [r, r*ell]), store.CRtV'];
        
        % Compute the relative residual norm
        % perp_egrad = svd_fact_approx(perp_egrad, max(3*r_increase, 20));
        perp_egrad = svd_fact(perp_egrad);
        % err_frob_exact = stable_norm_fact(residual) / problem.rhs_norm;
        err_frob = hutch_normest(residual, 10) / problem.rhs_norm;
        if err_frob > tol_resnorm % If not converged, increase the rank
            X0 = increase_rank(X, perp_egrad, r_increase);
            r_curr = r_curr + r_increase;
        end
        time_update = toc(timetic);
    end
    
    % Update the info and number of iterations
    info_tot = update_info(info_tot, info, time_update);
    iter = length(info_tot);

    % Reset slope history
    slope_history = zeros(1, 1000);
end

    function options = mergeDefaultOptions(options)
        % mergeDefaultOptions -  Merges the option structure given by the user with
        %   the default option structure. User-defined options have precedence
        %   over the default ones.
        %
        % Inputs:
        %   - options: Options structure with field options.optionname, where
        %       optionname is one of the following and the default value is indicated
        %       between parentheses:
        %           - r_start (1): starting rank r_0
        %           - r_increase (1): rank increase r_up
        %           - opts_solver: option structure for the fixed-rank solver
        %           - tol_resnorm (1e-6): tolerance on the residual error norm
        %           - maxiter_tot (10000): maximum total number of iterations
        %           - stopfun_options: options structure for heuristic
        %             criterions for stopping fixed-rank Riemannian
        %             optimization. The fields are:
        %               - miniter (5): minimum number of iterations before
        %                 triggering heuristic stopping criterions
        %               - patience (2): patience for heuristic criterion for
        %                 heuristic criterion based on sufficient decrease of
        %                 the gradient's norm
        %               - gradnorm_decrease (2): reduction factor for heuristic
        %                 criterion based on sufficient decrease of the
        %                 gradient's norm
        %               - backward_window (3): window to compute the slopes for
        %                 heuristic criterion based on slope of estimated residual
        %               - resnorm_slopefactor (0.5): tolerance factor for
        %                 heuristic criterion based on slope of estimated residual
        %               - epsilon_sigma (1e-10): relative tolerance for
        %                 truncating the rank of the iterate
        
        default_options = struct();
        default_options.r_start = 1;
        default_options.r_increase = 1;
        default_options.opts_solver = struct();
        default_options.tol_resnorm = 1e-6;
        default_options.maxiter_tot = 10000;
        stopfun_opts = struct();
        stopfun_opts.miniter = 5;
        stopfun_opts.patience = 2;
        stopfun_opts.backward_window = 3;
        stopfun_opts.resnorm_slopefactor = 0.5;
        stopfun_opts.gradnorm_decrease = 2;
        stopfun_opts.epsilon_sigma = default_options.tol_resnorm * 1e-4;
        default_options.stopfun_options = stopfun_opts;
        if isfield(options, "stopfun_options")
            options.stopfun_options = mergeOptions(default_options.stopfun_options, options.stopfun_options);
        end
        % Merge options with defaults
        options = mergeOptions(default_options, options);
    end

    function update_opts_solver()
        % update_opts_solver - Updates the options structure of the fixed-rank
        % Riemannian solver
        
        opts_solver.maxiter = min(opts_solver.maxiter, maxiter_tot - iter);
        if ~isempty(stopfun_options)
            opts_solver.stopfun = @fixed_rank_stopfun;
        end
        
        
        function stopnow = fixed_rank_stopfun(problem, x, info, last)
            % fixed_rank_stopfun - Checks the stopping criterions for
            % fixed-rank Riemannian optimization
            
            stopnow = false;
            miniter = stopfun_options.miniter;
            % Stopping criterion based on rank-deficient matrices
            if last >= miniter && trunc_rank(diag(x.S)) < r_curr
                stopnow = true;
            end
            
            % Stopping criterion based on sufficient decrease of gradient norm
            if ~stopnow
                patience = stopfun_options.patience;
                gradnorm_decrease = stopfun_options.gradnorm_decrease;
                if last >= miniter
                    suff_decrease = info(last).gradnorm < ...
                        (1/gradnorm_decrease) * max([info(last-patience:last-1).gradnorm]);
                    stopnow = ~suff_decrease;
                end
            end
            flag_stopfun = stopnow;
            
            % Stopping criterion based on slope
            if ~stopnow
                resnorm_slopefactor = stopfun_options.resnorm_slopefactor;
                backward_window = stopfun_options.backward_window;
                if last > backward_window
                    mean_slope = mean(slope_history(backward_window:last-1));
                else
                    mean_slope = 0;
                end
                if last >= backward_window
                    coeff = polyfit(1:backward_window, log([info(last-backward_window+1:last).res_norm_rel_est]), 1);
                    slope = coeff(1);
                    slope_history(last) = slope;
                    if last >= miniter
                        slopemin = min(slope_history(last-patience+1:last));
                        suff_slope = slopemin < resnorm_slopefactor * mean_slope;
                        stopnow = ~suff_slope;
                    end
                else
                    slope_history(last) = 0;
                end
            end
            flag_stopfun = stopnow;
        end
    end

    function k = trunc_rank(s)
        % trunc_rank - Computes the truncation rank based on a relative
        % criterion on the singular values
        %
        % Syntax:
        %   [k] = trunc_rank(s)
        %
        %
        % Inputs:
        %   - s: Vector of the singular values
        %
        % Outputs:
        %   - k: Truncation rank
        %
        
        sq_abs_err = cumsum([s(2:end).^2; 0], "reverse");
        sq_sum = sum(s.^2);
        sq_err = sq_abs_err / sq_sum;
        k = find(sq_err <= epsilon_sigma^2, 1);
    end

    function X0 = increase_rank(X, perp_egrad, r_up)
        % increase_rank - Computes a warm start for Riemannian optimization at
        % rank r + r_up, starting from the iterate X (with rank r) and
        % employing a normal correction step.
        %
        % Syntax:
        %   [X0] = increase_rank(X, perp_egrad, r_up)
        %
        % Inputs:
        %   - X         : Current iterate, with rank r
        %   - perp_egrad: Normal component of the Eucliden gradient
        %   - r_up      : Rank update
        %
        % Outputs:
        %   - X0: Starting point for optimization at rank r + r_up
        %
        
        Y = perp_egrad; Y.U = -Y.U;
        Y.U = Y.U(:, 1:r_up);
        Y.V = Y.V(:, 1:r_up);
        Y.S = Y.S(1:r_up, 1:r_up);
        if metric_EXD
            Y.EU = problem.M.E * Y.U;
            Y.DV = problem.M.D * Y.V;
        end
        
        % Compute the step size
        Ay = sylv_op_lr(problem.A, problem.B, Y);
        alpha = sum(diag(Y.S).^2) / mat_inner(Ay.U' * (Y.U * Y.S), Ay.V' * Y.V);
        assert(alpha > 0)
        
        % Update the rank
        X0.U = [X.U, Y.U]; X0.V = [X.V, Y.V]; X0.S = blkdiag(X.S, alpha * Y.S);
        if metric_EXD
            X0.EU = [X.EU, Y.EU]; X0.DV = [X.DV, Y.DV];
        end
        % Sort the SVD
        [~, idx] = sort(diag(X0.S), 'descend');
        X0.U = X0.U(:, idx); X0.V = X0.V(:, idx); X0.S = X0.S(idx, idx);
        if metric_EXD
            X0.EU = X0.EU(:, idx); X0.DV = X0.DV(:, idx);
        end
    end
    
    function Y = svd_fact(Y)
        % svd_fact - Computes the thin SVD factorization of Y, 
        % possibly with the nonstandard inner products induced by E and D.
        %
        % Syntax:
        %   [Y] = svd_fact(Y)
        %
        % Inputs:
        %   - Y: Structure with fields U and V, representing Y.U*Y.V'
        %
        % Outputs:
        %   - Y: Structure with fields U, S, V which are the thin SVD
        %        decomposition of the same matrix represented by the structure
        %        in input
        %
        
        if metric_EXD % Approximate E,D-SVD
            [Qu, Ru] = problem.M.qr_E(Y.U);
            [Qv, Rv] = problem.M.qr_D(Y.V);
        else % Approximate SVD
            [Qu, Ru] = qr(Y.U, 0);
            [Qv, Rv] = qr(Y.V, 0);
        end
        [U, S, V] = svd(Ru*Rv');
        Y.U = Qu * U;
        Y.V = Qv * V;
        Y.S = S;
    end


    function Y_approx = svd_fact_approx(Y, rsvd)
        % svd_fact_approx - Computes the thin approximate SVD factorization of Y, 
        % possibly with the nonstandard inner products induced by E and D.
        %
        % Syntax:
        %   [Y] = svd_fact_approx(Y)
        %
        % Inputs:
        %   - Y: Structure with fields U and V, representing Y.U*Y.V'
        %
        % Outputs:
        %   - Y: Structure with fields U, S, V which are the thin SVD
        %        decomposition of the same matrix represented by the structure
        %        in input
        %
        
        if metric_EXD % Approximate E,D-SVD
            [Qu, Ru] = qb(full(problem.M.cholE * Y.U), rsvd);
            Qu = problem.M.cholE \ Qu;

            [Qv, Rv] = qb(full(problem.M.cholD * Y.V), rsvd);
            Qv = problem.M.cholD \ Qv;
        else % Approximate SVD
            [Qu, Ru] = qb(Y.U, rsvd);
            [Qv, Rv] = qb(Y.V, rsvd);
        end
        [U, S, V] = svd(Ru*Rv');
        Y_approx.U = Qu * U;
        Y_approx.V = Qv * V;
        Y_approx.S = S;
    end

    function info_tot = update_info(info_tot, info, time_update)
        % update_info - Merges the info structures.
        
        % Add field with slope history
        aux = num2cell(slope_history(1:length(info))); 
        [info.slope_history] = aux{:};

        if isempty(info_tot)
            info(end).time = info(end).time + time_update;
            info_tot = info;
        else
            % Compute the total time
            aux = num2cell([info.time] + info_tot(end).time);
            [info.time] = aux{:};
            info(end).time = info(end).time + time_update;
            
            % Compute the total number of iterations
            aux = num2cell([info.iter] + info_tot(end).iter + 1);
            [info.iter] = aux{:};
            
            % Merge the structs
            info_tot = [info_tot, info];
        end
    end
end

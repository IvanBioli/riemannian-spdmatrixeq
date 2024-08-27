function [X, info] = pcg_trunc_lr(A, M_inv, C, T, options)
% pcg_trunc_lr - Preconditioned Conjugate Gradient with truncation
%
% Syntax:
%   [X, info] = pcg_trunc_lr(A, M_inv, C, T, options)
%   [X, info] = pcg_trunc_lr(A, M_inv, C, T)
%
% Description:
%   Solve the linear system A(X) = F via preconditioned conjugate gradient
%   with truncation, keeping the iterates in factored form. All linear
%   operators act on and return iterates in factored form.
%
% Inputs:
%   - A      : Linear operator
%   - M_inv  : Inverse of the preconditioner, i.e. M_inv ~ A^(-1)
%   - C      : Right-hand-side. See documentation of initialize
%   - T      : Truncation operator
%   - options: Option structure
%
% Outputs:
%   - X   : Approximate solution in factored form X.
%   - info: struct-array containing information about the iterations
%

options = mergeDefaultOptions(options);
if options.verbose
    fprintf("Running pCG with truncation\n")
    fprintf(' iter\t rank\t      rel. res\n');
end
% Get parameters from options
maxiter = options.maxiter;
tol = options.tol;
maxrank = options.maxrank;
truncate_R = options.truncate_R; opts_trunc_R = options.opts_trunc_R;
truncate_Q = options.truncate_Q; opts_trunc_Q = options.opts_trunc_Q;
truncate_P = options.truncate_P; opts_trunc_P = options.opts_trunc_P;
opts_trunc_X = options.opts_trunc_X;

% Initialization
timetic = tic(); % Time initialization
[X, rhs] = initialize(C, options);
iter = 0; rank_X = get_rank(X);
norm_rhs = stable_norm_fact(rhs);
R = rhs;
res = stable_norm_fact(R);
rel_res = res / norm_rhs;
% Save stats in a struct array info, and preallocate.
stats = savestats();
info(1) = stats;
info(min(10000, options.maxiter+1)).iter = [];
% Start iteration and timing
timetic = tic(); % Time the current iteration
Z = M_inv(R);
P = Z;
Q = A(P);
xi = inner_fact(P, Q);
update_trunc_tols();

% While loop
while iter < maxiter
    iter = iter+1;
    omega = inner_fact(R, P) / xi;
    X = lincomb_fact(X, omega, P);
    X = T(X, opts_trunc_X);
    rank_X = get_rank(X);
    
    stopnow = stopping_criterion(); if stopnow; break; end % Stop if we exceed the maximum rank
    
    AX = A(X);
    R = lincomb_fact(rhs, -1, AX);
    if truncate_R
        [R, res] = T(R, opts_trunc_R);
    else
        res = stable_norm_fact(R);
    end
    rel_res = res / norm_rhs;
    
    stopnow = stopping_criterion(); if stopnow; break; end % Stop if the residual is sufficiently low
    
    % Save the iteration infos
    stats = savestats();
    info(iter+1) = stats;
    timetic = tic(); % Time the current iteration
    
    Z = M_inv(R);
    beta = -inner_fact(Z, Q) / xi;
    P = lincomb_fact(Z, beta, P);
    if truncate_P
        P = T(P, opts_trunc_P);
    end
    Q = A(P);
    if truncate_Q
        Q = T(Q, opts_trunc_Q);
    end
    xi = inner_fact(P, Q);
end
info = info(1:iter+1);
if options.verbose
    if rel_res < tol
        fprintf("Convergence at iteration %d (relative error below" + ...
            " tolerance). options.tol = %e\nFinal relative error: %e\n", iter, tol, rel_res);
    elseif iter == maxiter
        fprintf("Reached maximum number of iterations; options.maxiter = %d" + ...
            "\nFinal relative error: %.2e\n", maxiter, rel_res);
    elseif rank_X > maxrank
        fprintf("Reached maximum rank; options.maxrank = %d" + ...
            "\nFinal relative error: %e\n\n", maxrank, rel_res);
    else
        fprintf("Detected stagnation" + ...
            "\nFinal relative error: %.2e\n", maxiter, rel_res);
    end
    fprintf('Total time is %f [s]\n', info(end).time);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AUXILIARY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function update_trunc_tols()
        % update_trunc_tols - Defines the absolute and relative the
        % truncation tolerances for R, P, Q
        
        if ~isempty(opts_trunc_R.method) && ismember("absF_fact", opts_trunc_R.method)
            i = find(opts_trunc_R.method == "absF_fact");
            opts_trunc_R.method(i) = "absF";
            opts_trunc_R.tol(i) = opts_trunc_R.tol(i) * norm_rhs;
        end
        
        if ~isempty(opts_trunc_P.method) &&ismember("absF_fact", opts_trunc_P.method)
            i = find(opts_trunc_P.method == "absF_fact");
            opts_trunc_P.method(i) = "absF";
            opts_trunc_P.tol(i) = opts_trunc_P.tol(i) * stable_norm_fact(P);
        end
        
        if ~isempty(opts_trunc_Q.method) && ismember("absF_fact", opts_trunc_Q.method)
            i = find(opts_trunc_Q.method == "absF_fact");
            opts_trunc_Q.method(i) = "absF";
            opts_trunc_Q.tol(i) = opts_trunc_Q.tol(i) * stable_norm_fact(Q);
        end
        
    end

    function stopnow = stopping_criterion()
        % stopping_criterion - Returns whether any of the prescribed
        % stopping criterion is met
        
        stopnow = false;
        % Stopping criterion based on relative tolerance or maxrank
        if  (rel_res < tol) || (rank_X >= maxrank)
            stopnow = true;
        end
        
        % Stopping criterion based on number of iterations
        if  iter >= maxiter
            stopnow = true;
        end
        
        % Stopping criterion to avoid stagnation
        if options.detect_stagnation && iter > 3
            if rel_res >  (max([info(iter-3:iter-1).res_norm_rel]) + tol * 1e-2)
                stopnow = true;
            end
        end
        
        % Log info before stopping
        if stopnow
            stats = savestats();
            info(iter+1) = stats;
        end
        
    end

    function stats = savestats()
        % savestats - Routine in charge of collecting the current iteration
        % stats.
        
        stats.iter = iter;
        stats.res_norm = res;
        stats.res_norm_rel = rel_res;
        if iter == 0
            stats.rank = 0;
            stats.time = toc(timetic);
        else
            stats.rank = rank_X;
            stats.time = info(iter).time + toc(timetic);
        end
        % Print logs
        if options.verbose
            fprintf('%5d\t%5d\t%.8e\n', iter, rank_X, rel_res);
        end
    end
end

function X = lincomb_fact(X1, alpha, X2)
% lincomb_fact - Computes the linear combination X = X1 + alpha * X2 in
% factored form.

if isfield(X1, 'V') && ~isfield(X1, 'U')
    X.V = [X1.V, X2.V];
    X.D = blkdiag(X1.D, alpha * X2.D);
elseif isfield(X1, 'U')
    X.U = [X1.U, alpha*X2.U];
    X.V = [X1.V, X2.V];
elseif ~isstruct(X1)
    X = X1 + alpha * X2;
else
    error("Unrecognized struct for X")
end
end

function inner = inner_fact(X1, X2)
% inner_fact - Compute the Frobenius inner product of X1 and X2

if isfield(X1, 'V') && ~isfield(X1, 'U')
    aux = X2.V' * X1.V;
    fact_1 = X2.D * aux;
    fact_2 = aux * X1.D;
    inner = fact_1(:)'*fact_2(:);
elseif isfield(X1, 'U')
    fact_1 = X2.U' * X1.U;
    fact_2 = X2.V' * X1.V;
    inner = fact_1(:)'*fact_2(:);
elseif ~isstruct(X1)
    inner = X1(:)'*X2(:);
else
    error("Unrecognized struct for X")
end
end

function r = get_rank(X)
% get_rank - Return the numerical rank of the current iterate

if isfield(X, 'V')
    r = size(X.V, 2);
elseif ~isstruct(X)
    r = rank(X);
else
    error("Unrecognized struct for X")
end
end

function [X, rhs] = initialize(C, options)
% initialize - Defines the structure of the initial iterate, based on the 
% format of C, i.e. the right-hand-side.
%
% Inputs:
%   - C      : Structure of the right-hand-side. It can be:
%       - structure with fields L and R. In this case the RHS is F =
%         C.L*C.R', and iterates are kept in the low-rank format X.U * X.V'
%       - matrix. If options.lowrank, the RHS is F = C*C' and iterates are 
%         kept in the low-rank format X.V * X.D * X.V'. Otherwise, the 
%         iterates are stored as full matrices and F = C. 
%   - options: Options structure
%
% Outputs:
%   - X  : initial iterate structure
%   - rhs: right-hand side structure
%

if ~isstruct(C) && options.lowrank
    rhs.V = C; rhs.D = speye(size(C, 2));
    X.V = []; X.D = []; % Initializing X to 0
elseif ~isstruct(C) && ~options.lowrank
    rhs = C;
    X = zeros(size(C));
else
    rhs.U = C.L; rhs.V = C.R;
    X.U = []; X.V = [];
end
end

function options = mergeDefaultOptions(options)
% mergeDefaultOptions - Merges the option structure given by the user with 
%   the default option structure. User-defined options have precedence 
%   over the default ones. 
%
% Inputs:
%   - options: Options structure with field options.optionname, where
%       optionname is one of the following and the default value is indicated
%       between parentheses:
%           - maxiter (100): maximum number of iterations
%           - tol (1e-6): relative tolerance on the residual
%           - lowrank (true): whether iterates are in low-rank format
%           - maxrank (Inf): maximum rank of intermediate iterates
%           - truncate_R (true): whether to truncate the R
%           - truncate_Q (false): whether to truncate the Q
%           - truncate_P (true): whether to truncate the P
%           - verbose (true): whether to print to screen convergence logs
%           - detect_stagnation (true): whether to stop iterations when
%             stagnation of the residual is detected
%           - opts_trunc_R: option structure for truncation of R
%           - opts_trunc_P: option structure for truncation of P
%           - opts_trunc_Q: option structure for truncation of Q
%           - opts_trunc_X: option structure for truncation of X
%
% Outputs:
%   - options: Options structure containing all the fields above, with
%     default values in fields that have not been provided by the user
%

% Define default options
default_options.maxiter = 100;
default_options.tol = 1e-6;
default_options.lowrank = true;
default_options.maxrank = Inf;
default_options.truncate_R = true;
default_options.truncate_Q = false;
default_options.truncate_P = true;
default_options.verbose = true;
default_options.detect_stagnation = false;
default_options.opts_trunc_R = struct("method", [], "tol", []);
default_options.opts_trunc_P = struct("method", [], "tol", []);
default_options.opts_trunc_Q = struct("method", [], "tol", []);
default_options.opts_trunc_X = struct("method", [], "tol", []);


% Merge data with defaults
options = mergeOptions(default_options, options);
end
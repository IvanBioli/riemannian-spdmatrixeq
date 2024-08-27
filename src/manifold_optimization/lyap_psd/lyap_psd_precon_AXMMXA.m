function [g, store] = lyap_psd_precon_AXMMXA(X, eta, A, store, options, M)
% lyap_psd_precon_AXMMXA - Applies to eta the
% Riemannian preconditioner derived from the generalized Lyapunov operator
% P(X) = AXM + MXA, i.e. solves the projected generalized Lyapunov
% equation Proj_X (A * g * M + M * g * A) = eta, for g in Tang_X(Snr++). If
% M is specified, it assumes that we are employing the Riemannian
% optimization tools for M_r as an embedded sumbmanifol of Rmn with the
% metric induced by X -> MXM.
%
% Syntax:
%   [g, store] = lyap_psd_precon_AXMMXA(X, eta, A, store, options, M)
%
% Inputs:
%   - X      : Current iterate
%   - eta    : Vector of the tangent space at X to which the preconditioner
%              is applied
%   - A      : Coefficient matrix of the preconditioning operator
%   - store  : Manopt's StoreDB struct for the current iterate
%   - options: Options structure.  See mergeDefaultOptions
%   - M      : Coefficient matrix of the preconditioning operator. If not
%              specified it is assumed to be the identity matrix
%
% Outputs:
%   - g    : The result of the application of the preconditioner to eta
%   - store: Updated store struct
%

if nargin < 5
    options = struct();
end
options = mergeDefaultOptions(options);
is_M = (nargin >= 6);

% Computations independent of eta
[n, r] = size(X.V);
if is_M
    store = prepare(X, A, store, options, M);
end
if ~is_M
    store = prepare(X, A, store, options);
    M = speye(n);
end
MV = store.precond.MV_newbasis;

% Change of basis
Q = store.precond.Q;
eta.Vp = eta.Vp * Q;  eta.M = Q' * eta.M * Q;
if is_M
    eta.MVp = eta.MVp * Q;
else
    eta.MVp = eta.Vp;
end

% Compute W
W_step1 = zeros(n, 1, r);
for i = 1:r % FIXME: Embarassigly parallel
    W_step1(:, :, i) = (A + store.precond.lambdas(i) * M) \ eta.MVp(:, i);
end
if ~options.small_memory
    W = W_step1 - pagemtimes(store.precond.S_step1, ...
        pagemldivide(store.precond.Sneg, pagemtimes(MV', W_step1)));
    Wv = squeeze(W);
else
    Wv = squeeze(W_step1);
    Lambda = reshape(store.precond.lambdas, [1, 1, r]) .* repmat(eye(r), [1, 1, r]);
    for i = 1:r
        SFv_step1_i = (A + store.precond.lambdas(i) * M) \ store.precond.MvPerpAV;
        S_step1_i = SFv_step1_i(:, 1:r);
        Sneg_i = symmetrize(MV' * S_step1_i);
        Wv(:, i) = Wv(:, i) - S_step1_i * (Sneg_i \ (MV' * Wv(:, i)));
        Fv_step1_i = SFv_step1_i(:, r+1:end);
        Fv_i = Fv_step1_i - S_step1_i * (Sneg_i \ (MV' * Fv_step1_i));
        Lambda(:, :, i) = Lambda(:, :, i) - store.precond.AV' * Fv_i;
    end
    Lambda = symmetrize(Lambda);
    store.precond.Lambda = Lambda;
    store = prepare(X, A, store, options, M); % Assemble F
end

% Assemble the RHS
AV = store.precond.AV; AVtWv = AV' * Wv;
RHS = eta.M - AVtWv - (AVtWv)'; RHS = RHS(:);

% Assemble F (possibly only as an operator) and solve the system
if options.assemble_F % Solve the system directly
    vec_Mxi = store.precond.F \ RHS;
    %fprintf("\t norm(F * vec_Mxi - RHS) / norm(RHS): %e\n", norm(F * vec_Mxi - RHS) / norm(RHS))
else % Solve the system using CG
    [vec_Mxi, flag, rel_res, iter] = pcg(store.precond.F_fun, RHS, ...
        options.pcg_tol, options.pcg_maxit); %rel_res and iter useful only for debugging
    if options.pcg_verbose && flag >= 1
        warning('pcg failed to converge')
    end
    %fprintf("\t norm(F * vec_Mxi - RHS) / norm(RHS): %e\n", norm(F_fun(vec_Mxi) - RHS) / norm(RHS))
end

% Obtain Mxi
Mxi = reshape(vec_Mxi, [r, r]);
% Make numerically symmetric and throw warning
tol = 1e-4;
Mxi_old = Mxi; Mxi = .5*(Mxi + Mxi');
err_M = norm(Mxi - Mxi_old, 'fro')/ norm(Mxi_old, 'fro');
if err_M > tol
    warning("Rel. error on tangent condition for g.M:%e\n", err_M)
end
% Obtain Vxi
if ~options.small_memory
    Vxi = Wv - ...
        squeeze(pagemtimes(store.precond.Fv, reshape(Mxi, [r, 1, r])));
else
    Vxi = Wv;
    Mxi_reshaped = reshape(Mxi, [r, r]);
    for i = 1:r
        SFv_step1_i = (A + store.precond.lambdas(i) * M) \ store.precond.MvPerpAV;
        S_step1_i = SFv_step1_i(:, 1:r);
        Sneg_i = symmetrize(MV' * S_step1_i);
        Fv_step1_i = SFv_step1_i(:, r+1:end);
        Fv_i = Fv_step1_i - S_step1_i * (Sneg_i \ (MV' * Fv_step1_i));
        Vxi(:, i) = Vxi(:, i) - Fv_i * Mxi_reshaped(:, i);
    end
end

% Change back the basis
g.Vp = Vxi * Q';
g.M = Q * Mxi * Q';
V = X.V;
if is_M
    MV = X.MV;
else
    MV = V;
end

% Ensure that g.M and g.Vp numerically verify the tangent constraints
g_orth.M = 0.5 * (g.M + g.M');
g_orth.Vp = g.Vp - V*(MV'*g.Vp);
% Check and throw warning
err_Vp = norm(g.Vp - g_orth.Vp, 'fro')/ norm(g.Vp, 'fro');
if err_Vp > tol
    warning("Rel. error on tangent condition for g.Vp:%e\n", err_Vp)
end
% Set g to the numerically more stable version
g = g_orth;
end

function store = prepare(X, A, store, options, M)
% prepare - Performs the computations that are independent of the tangent
% vector eta, anche caches results that can be subsequently reused in the
% store structure

% Precomputations
if ~isfield(store, 'precond')
    store.precond = struct();
    % Definition of relevant quantities
    V = X.V;
    [n, r] = size(V);
    is_M = (nargin > 4);
    if ~is_M
        M = speye(n);
    end
end

% Change of basis
if ~isfield(store.precond, 'Q')
    AV = A * V; VtAV = V' * AV; VtAV = symmetrize(VtAV);
    % Change of basis
    [Q, lambdas] = eig(VtAV); lambdas = diag(lambdas);
    store.precond.Q = Q; store.precond.lambdas = lambdas;
end
if ~isfield(store.precond, 'MV_newbasis')
    V = V * Q;
    if is_M
        MV = X.MV * Q;
    else
        MV = V;
    end
    store.precond.MV_newbasis = MV;
end

% Compute Sv, Fv
if ~isfield(store.precond, 'AV')
    AV = A * V; PerpAV = AV - MV .* lambdas'; MVPerpAV = [MV, PerpAV];
    store.precond.AV = AV;
    store.precond.MvPerpAV = MVPerpAV;
end
if ~isfield(store.precond, 'Sneg') && ~options.small_memory
    % Solve the shifted linear systems
    SFv_step1 = zeros(n, 2*r, r);
    for i = 1:r % FIXME: Embarassigly parallel
        SFv_step1(:, :, i) = (A + lambdas(i) * M) \ MVPerpAV;
    end
    
    % Compute S
    S_step1 = SFv_step1(:, 1:r, :);
    Sneg = pagemtimes(MV', S_step1); Sneg = symmetrize(Sneg);
    store.precond.S_step1 = S_step1; store.precond.Sneg = Sneg;
    % Compute Fv
    Fv_step1 = SFv_step1(:, r+1:end, :);
    Fv = Fv_step1 - pagemtimes(S_step1, pagemldivide(Sneg, pagemtimes(MV', Fv_step1)));
    store.precond.Fv = Fv;
end

% Assemble F (possibly only as an operator)
if ~options.small_memory && ~isfield(store.precond, 'Lambda')
    Lambda = reshape(lambdas, [1, 1, r]) .* repmat(eye(r), [1, 1, r]) ...
        - pagemtimes(AV', Fv);
    Lambda = symmetrize(Lambda); store.precond.Lambda = Lambda;
end
if ~isfield(store.precond, 'F') && ~isfield(store.precond, 'F_fun') && isfield(store.precond, 'Lambda')
    Lambda = store.precond.Lambda;
    Lambda_cell = cellfun(@sparse,  num2cell(Lambda,[1,2]), 'uni',0);
    if options.assemble_F % Assemble F
        r = size(X.V,2); Pi_idx = TvecMat_indices(r, r);
        F = blkdiag(Lambda_cell{:}); F = F + F(Pi_idx, Pi_idx);
        F = symmetrize(F);
        if nnz(F) / numel(F) > 0.25
            F = full(F);
        end
        store.precond.F = decomposition(F, 'chol');
    else % Assemble F as an operator
        F_fun =@(X) F_matvec(X, Lambda, r);
        store.precond.F_fun = F_fun;
    end
end
end

function FX = F_matvec(X, Lambda, r)
% F_matvec - Applies to the matrix X the linear operator that constitutes
% the LHS of the final system to be solved to compute Mxi (see report).
% Exploits the matrix-matrix formulation and not the Kronecker one.

X = reshape(X, [r, r]);
FX = squeeze(pagemtimes(Lambda, reshape(X, [r, 1, r]))) + ...
    squeeze(pagemtimes(Lambda, reshape(X', [r, 1, r])))';
FX = FX(:);
end

function opts = mergeDefaultOptions(opts)
% mergeDefaultOptions - Merges the option structure given by the user with
%   the default option structure. User-defined options have precedence
%   over the default ones.
%
% Inputs:
%   - opts: Options structure with field opts.optionname, where optionname
%       is one of the following and the default value is indicated
%       between parentheses:
%           - assemble_F (true): whether to assemble explicitly the matrix
%             F of the final system, or to use the corresponding linear
%             operator and employ pCG
%           - pcg_tol (1e-6): (used only if assemble_F == false) Tolerance
%             for pCG
%           - pcg_maxit (250): (used only if assemble_F == false) Maximum
%             number of iterations for pCG
%           - pcg_verbose (true): whether to print if pCG converged
%           - small_memory (false): whether the O(r^2 n) memory requirement
%             is too much and caching should be avoided

default_opts.assemble_F = true;
default_opts.pcg_tol = 1e-6;
default_opts.pcg_maxit = 250;
default_opts.pcg_verbose = true;
default_opts.small_memory = false;
opts = mergeOptions(default_opts, opts);
end
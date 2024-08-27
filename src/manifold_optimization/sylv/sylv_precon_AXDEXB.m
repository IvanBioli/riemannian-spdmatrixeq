function [g, store] = sylv_precon_AXDEXB(X, eta, A, B, store, options, D, E)
% sylv_precon_AXDEXB - Applies to eta the Riemannian
% preconditioner derived from the generalized Sylvester operator
% P(X) = AXD + EXB, i.e. solves the projected generalized Sylvester
% equation Proj_X (A * g * D + E * g * B) = eta, for g in Tang_X(M_r). If
% D,E are specified, it assumes that we are employing the Riemannian
% optimization tools for M_r as an embedded sumbmanifol of Rmn with the
% metric induced by X -> EXD.
%
% Syntax:
%   [g, store] = sylv_precon_AXDEXB(X, eta, A, B, store, options, D, E)
%
% Inputs:
%   - X      : Current iterate
%   - eta    : Vector of the tangent space at X to which the preconditioner
%              is applied
%   - A, B   : Coefficient matrices of the preconditioning operator
%   - store  : Manopt's StoreDB struct for the current iterate
%   - options: Options structure.  See mergeDefaultOptions
%   - D      : Coefficient matrix of the preconditioning operator. If not
%              specified it is assumed to be the identity matrix
%   - E      : Coefficient matrix of the preconditioning operator. If not
%              specified it is assumed to be the identity matrix
%
% Outputs:
%   - g    : The result of the application of the preconditioner to eta
%   - store: Updated store struct
%

if nargin < 6
    options = struct();
end
options = mergeDefaultOptions(options);
is_ED = nargin >= 7;

% Computations independent of eta
[m, r] = size(X.U); [n, ~] = size(X.V);
if is_ED
    store = prepare(X, A, B, store, options, D, E);
else
    store = prepare(X, A, B, store, options);
    E = speye(m);
    D = speye(n);
end
EU = store.precond.EU_newbasis;
DV = store.precond.DV_newbasis;

% Change of basis
Qu = store.precond.Qu; Qv = store.precond.Qv;
eta.Up = eta.Up*Qu; eta.Vp = eta.Vp*Qv; eta.M = Qv' * eta.M * Qu;
if is_ED
    eta.EUp = eta.EUp * Qu; eta.DVp = eta.DVp * Qv;
else
    eta.EUp = eta.Up; eta.DVp = eta.Vp;
end

% Compute Wu and Wv
Wu_step1 = zeros(m, 1, r); Wv_step1 = zeros(n, 1, r);
for i = 1:r % FIXME: Embarassigly parallel
    if options.cache_decomp
        Wu_step1(:, :, i) = store.precond.AlambdaE{i} \ eta.EUp(:, i);
        Wv_step1(:, :, i) = store.precond.BkappaD{i} \ eta.DVp(:, i);
    else
        Wu_step1(:, :, i) = (A + store.precond.lambdas(i) * E) \ eta.EUp(:, i);
        Wv_step1(:, :, i) = (B + store.precond.kappas(i) * D) \ eta.DVp(:, i);
    end
end
Wu = Wu_step1 - pagemtimes(store.precond.Su_step1, ...
    pagemldivide(store.precond.Suneg, pagemtimes(EU', Wu_step1)));
Wv = Wv_step1 - pagemtimes(store.precond.Sv_step1, ...
    pagemldivide(store.precond.Svneg, pagemtimes(DV', Wv_step1)));
Wu = squeeze(Wu); Wv = squeeze(Wv);

% Assemble the RHS
RHS = eta.M - store.precond.AU' * Wu - (store.precond.BV' * Wv)';
RHS = RHS(:);

% Solve the system
if options.assemble_F % Solve the system using a direct solver
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

% Obtain Uxi and Vxi
Mxi = reshape(vec_Mxi, [r, r]);
Uxi = Wu - squeeze(pagemtimes(store.precond.Fu, reshape(Mxi, [r, 1, r])));
Vxi = Wv - squeeze(pagemtimes(store.precond.Fv, reshape(Mxi', [r, 1, r])));

% Change back the basis
g.Up = Uxi * Qu';
g.Vp = Vxi * Qv';
g.M = Qv * Mxi * Qu';

% Ensure that g.M and g.Vp numerically verify the tangent constraints
V = X.V; U = X.U;
if is_ED
    EU = X.EU; DV = X.DV;
else
    EU = U; DV = V;
end
g_orth.M = g.M;
g_orth.Vp = g.Vp - V*(DV'*g.Vp);
g_orth.Up = g.Up - U*(EU'*g.Up);
% Check and throw warning
tol = 1e-6;
err_Vp = norm(g.Vp - g_orth.Vp, 'fro')/ norm(g_orth.Vp, 'fro');
if err_Vp > tol
    warning("Rel. error on tangent condition for g.Vp:%e\n", err_Vp)
end
err_Up = norm(g.Up - g_orth.Up, 'fro')/ norm(g_orth.Up, 'fro');
if err_Up > tol
    warning("Rel. error on tangent condition for g.Up:%e\n", err_Up)
end
% Set g to the numerically more stable version
g = g_orth;
end

function store = prepare(X, A, B, store, options, D, E)
% prepare - Performs the computations that are independent of the tangent
% vector eta, anche caches results that can be subsequently reused in the
% store structure

% Precomputations
if ~isfield(store, 'precond')
    store.precond = struct();
    % Definition of relevant quantities
    U = X.U; V = X.V;
    [m, r] = size(U); [n, ~] = size(V);
    is_ED = nargin > 5;
    if ~is_ED
        E = speye(m);
        D = speye(n);
    end
    Pi_idx = TvecMat_indices(r, r);
end

% Change of basis
if ~isfield(store.precond, 'Qu')
    BV = B * V; VtBV = V' * BV; VtBV = symmetrize(VtBV);
    [Qu, lambdas] = eig(VtBV); lambdas = diag(lambdas);
    store.precond.Qu = Qu; store.precond.lambdas = lambdas;
    
    AU = A * U; UtAU = U' * AU; UtAU = symmetrize(UtAU);
    [Qv, kappas] = eig(UtAU); kappas = diag(kappas);
    store.precond.Qv = Qv; store.precond.kappas = kappas;
end

% Change of basis
if ~isfield(store.precond, 'EU_newbasis')
    V = V * Qu;
    U = U * Qv;
    if is_ED
        EU = X.EU * Qv; DV = X.DV * Qu;
    else
        EU = U; DV = V;
    end
    store.precond.EU_newbasis = EU;
    store.precond.DV_newbasis = DV;
end

% Compute Su, Sv, Fu, Fv
if ~isfield(store.precond, 'Suneg')
    AU = A * U; PerpAU = AU - EU .* kappas'; EUPerpAU = [EU, PerpAU];
    BV = B * V; PerpBV = BV - DV .* lambdas'; DVPerpBV = [DV, PerpBV];
    store.precond.AU = AU; store.precond.BV = BV;
    % Solve the shifted linear systems
    SFu_step1 = zeros(m, 2*r, r); SFv_step1 = zeros(n, 2*r, r);
    if options.cache_decomp
        store.precond.AlambdaE = cell(r, 1);
        store.precond.BkappaD = cell(r, 1);
    end
    for i = 1:r % FIXME: Embarassigly parallel
        if options.cache_decomp
            store.precond.AlambdaE{i} = decomposition(A + lambdas(i) * E);
            SFu_step1(:, :, i) = store.precond.AlambdaE{i} \ EUPerpAU;
            store.precond.BkappaD{i} = decomposition(B + kappas(i) * D);
            SFv_step1(:, :, i) = store.precond.BkappaD{i} \ DVPerpBV;
        else
            SFu_step1(:, :, i) = (A + lambdas(i) * E) \ EUPerpAU;
            SFv_step1(:, :, i) = (B + kappas(i) * D) \ DVPerpBV;
        end
    end
    
    % Compute Su, Sv
    Su_step1 = SFu_step1(:, 1:r, :); Sv_step1 = SFv_step1(:, 1:r, :);
    Suneg = pagemtimes(EU', Su_step1); Svneg = pagemtimes(DV', Sv_step1);
    Suneg = symmetrize(Suneg); Svneg = symmetrize(Svneg);
    store.precond.Su_step1 = Su_step1; store.precond.Sv_step1 = Sv_step1;
    if options.cache_decomp
        Suneg = pagedecomposition(Suneg, 'chol');
        Svneg = pagedecomposition(Svneg, 'chol');
    end
    store.precond.Suneg = Suneg; store.precond.Svneg = Svneg;
    % Compute Fu and Fv
    Fu_step1 = SFu_step1(:, r+1:end, :); Fv_step1 = SFv_step1(:, r+1:end, :);
    Fu = Fu_step1 - pagemtimes(Su_step1, pagemldivide(Suneg, pagemtimes(EU', Fu_step1)));
    Fv = Fv_step1 - pagemtimes(Sv_step1, pagemldivide(Svneg, pagemtimes(DV', Fv_step1)));
    store.precond.Fu = Fu; store.precond.Fv = Fv;
end

% Assemble F (possibly only as an operator)
if ~isfield(store.precond, 'F') && ~isfield(store.precond, 'F_fun')
    Lambda = reshape(lambdas, [1, 1, r]) .* repmat(eye(r), [1, 1, r]) ...
        - pagemtimes(AU', Fu);
    Lambda = symmetrize(Lambda);
    Lambda_cell = cellfun(@sparse,  num2cell(Lambda,[1,2]), 'uni',0);
    Kappa = reshape(kappas, [1, 1, r]) .* repmat(eye(r), [1, 1, r]) ...
        - pagemtimes(BV', Fv);
    Kappa = symmetrize(Kappa);
    Kappa_cell = cellfun(@sparse,  num2cell(Kappa,[1,2]), 'uni',0);
    if options.assemble_F % Assemble F
        F = blkdiag(Kappa_cell{:}); F = F(Pi_idx, Pi_idx);
        F = F + blkdiag(Lambda_cell{:});
        F = symmetrize(F);
        store.precond.F = decomposition(F, 'chol');
    else % Assemble F as an operator
        F_fun =@(X) F_matvec(X, Kappa, Lambda, r);
        store.precond.F_fun = F_fun;
    end
end
end

function FX = F_matvec(X, Kappa, Lambda, r)
% F_matvec - Applies to the matrix X the linear operator that constitutes
% the LHS of the final system to be solved to compute Mxi (see report).
% Exploits the matrix-matrix formulation and not the Kronecker one.

X = reshape(X, [r, r]);
FX = squeeze(pagemtimes(Lambda, reshape(X, [r, 1, r]))) + ...
    squeeze(pagemtimes(Kappa, reshape(X', [r, 1, r])))';
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
%           - cache_decomp (false): whether to cache the decomposition of
%             the shifted matrices and of the Schur complements.
%           - assemble_F (true): whether to assemble explicitly the matrix
%             F of the final system, or to use the corresponding linear
%             operator and employ pCG
%           - pcg_tol (1e-6): (used only if assemble_F == false) Tolerance
%             for pCG
%           - pcg_maxit (250): (used only if assemble_F == false) Maximum
%             number of iterations for pCG
%           - pcg_verbose (true): whether to print if pCG converged

default_opts.cache_decomp = false;
if isfield(opts, 'cache_decomp' ) && opts.cache_decomp
    error("Caching decomposition proved not to be effective in practice")
end
default_opts.assemble_F = true;
default_opts.pcg_tol = 1e-6;
default_opts.pcg_maxit = 250;
default_opts.pcg_verbose = true;
opts = mergeOptions(default_opts, opts);
end
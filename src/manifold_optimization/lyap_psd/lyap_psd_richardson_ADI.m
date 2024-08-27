function [g, store] = lyap_psd_richardson_ADI(X, A, M, N, B, store, options)
% lyap_psd_richardson_ADI - Applies the Riemannian truncated
% preconditioned Richardson iteration with generalized-Lyapunov-like
% preconditioner. The inverse of the generalized Lyapunov operator is
% approximated by steps of fADI algorithm or of the bilinear ADI algorithm.
%
% Syntax:
%   [g, store] = lyap_psd_richardson_ADI(X, A, M, N, B, store, options)
%
% Inputs:
%   - X      : Current iterate
%   - A, M   : Coefficient matrices of the multiterm linear matrix equation
%   - N      : Cell array containing the other coefficient matrices of the
%              multiterm linear matrix equation
%   - B      : Factor of the righ-hand-side B*B'
%   - store  : Manopt's StoreDB struct for the current iterate
%   - options: Options structure. It is a structure with fields:
%                - fadi: structure with fields A, M, N, p. The matrices
%                  A, M and the cell array N define the preconditioning
%                  operator. If N is present, the preconditioner is
%                       P(X) = AXM + MXA - N{1} X N{1} - ... - N{k} X N{k}'
%                   and the bilinear ADI iteration is employed. Otherwise,
%                   the preconditioning operator is P(X) = AXM + MXA and
%                   the fADI iteration is employed. The vector p contains
%                   the ADI shift parameters. If M is not specified it is
%                   assumed to be the identity matrix.
%                - trunc_opts: structure with truncation options for
%                  truncating the Euclidean gradient before applying fADI.
%
% Outputs:
%   - g    : Preconditioned search direction
%   - store: Updated store struct
%

% Precomputations
store = lyap_psd_prepare_cost(X, A, M, N, B, store);
[n, r] = size(X.V);
ell = length(N);

% Assemble L(X) - C in factored form
g = struct();
g.V = [(store.AV + store.MV) / sqrt(2), ...
    (store.AV - store.MV) / sqrt(2), ...
    reshape(store.NV, [n, r*ell]), ...
    B];
g.D = blkdiag(kron(...
    spdiags([1;-1;-ones(ell, 1)], 0, ell+2, ell+2), ...
    X.D), ...
    -speye(size(B,2)));

% Truncate g
if isfield(options, "trunc_opts")
    g = truncation(g, options.trunc_opts);
end

% Apply ADI iteration, possibly bilinear
if ~isfield(options.fadi, 'options')
    options.fadi.options = struct();
end
if ~isfield(options.fadi, 'M')
    options.fadi.M = [];
end
if isfield(options.fadi, 'N')
    g_precond = adi_lyap(options.fadi.A, g, options.fadi.p, ...
        options.fadi.options, options.fadi.M, options.fadi.N);
else
    g_precond = adi_lyap(options.fadi.A, g, options.fadi.p, ...
        options.fadi.options, options.fadi.M);
end

% Project to the tangent
g = project(X, g_precond);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AUXILIARY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Zproj = project(X, Z)
% project - Projects the preconditioned Euclidean gradient to the tangent
% space

symZV = Z.V*(Z.D*(Z.V'*X.V));
VtsymZV = X.V'*symZV;
Zproj.M = VtsymZV;
Zproj.Vp = symZV  - X.V*VtsymZV;
end


function store = sylv_prepare_cost(X, A, B, C, store)
% sylv_prepare_cost - Performs the precomputations for cost
% and gradient, with the energy-norm based cost function.
%
% Syntax:
%   [store] = sylv_prepare_cost(X, A, B, C, store)
%
%
% Inputs:
%   - X    : Current iterate
%   - A, B : Cell arrays containing the coefficient matrices of the
%            multiterm linear matrix equation
%   - C    : Structure of the right-hand-side. It is a structure with
%            fields L and R, and the RHS is C.L * C.R'
%   - store: Manopt's StoreDB struct for the current iterate
%
% Outputs:
%   - store: Updated store structure, with cached precomputations
%

% Precomputations
U = X.U; [m, r] = size(U);
V = X.V; [n, ~] = size(V);
s = diag(X.S);
ell = length(A);

% Compute the products A * U, B * V, C_R' * V and C_L' * U
if ~isfield(store, 'AU')
    store.AU = zeros(m, r, ell);
    for i = 1:ell
        store.AU(:, :, i) = A{i} * U;
    end
end
if ~isfield(store, 'BV')
    store.BV = zeros(n, r, ell);
    for i = 1:ell
        store.BV(:, :, i) = B{i}*V;
    end
end
if ~isfield(store, 'CLtU')
    store.CLtU = C.L' * U;
end
if ~isfield(store, 'CRtV')
    store.CRtV = C.R' * V;
end

% Compute the products U' * A * U and V' * B * V
if ~isfield(store, 'UtAU')
    store.UtAU = pagemtimes(U', store.AU);
end
if ~isfield(store, 'UtAUS')
    store.UtAUS = bsxfun(@times, store.UtAU, s');
end
if ~isfield(store, 'UtAtUS')
    store.UtAtUS = bsxfun(@times, pagetranspose(store.UtAU), s');
end
if ~isfield(store, 'VtBV')
    store.VtBV = pagemtimes(V', store.BV);
end
if ~isfield(store, 'SVtBV')
    store.SVtBV = bsxfun(@times, s, store.VtBV);
end
if ~isfield(store, 'SVtBtV')
    store.SVtBtV = bsxfun(@times, s, pagetranspose(store.VtBV));
end
end
function store = lyap_psd_prepare_cost(X, A, M, N, B, store)
% lyap_psd_prepare_cost - Performs the precomputations for cost
% and gradient, with the energy-norm based cost function.
%
% Syntax:
%   [store] = lyap_psd_prepare_cost(X, A, M, N, B, store)
%
% Inputs:
%   - X    : Current iterate
%   - A, M : Coefficient matrices of the multiterm linear matrix equation
%   - N    : Cell array containing the other coefficient matrices of the
%              multiterm linear matrix equation
%   - B    : Factor of the righ-hand-side B*B'
%   - store: Manopt's StoreDB struct for the current iterate
%
% Outputs:
%   - store: Updated store structure, with cached precomputations
%


% Get the number of terms
ell = length(N);

% Get the fields of X
d = diag(X.D);
V = X.V;
[n, r] = size(V);

is_M = ~isempty(M);

% Compute the products A * V, M * V and N_i * V
if ~isfield(store, 'AV')
    store.AV = A * V;
end
if ~isfield(store, 'MV')
    if is_M
        store.MV = M * V;
    else
        store.MV = V;
    end
end
if ~isfield(store, 'NV')
    store.NV = zeros(n, r, ell);
    for i = 1:ell
        store.NV(:, :, i) = N{i} * V;
    end
end
if ~isfield(store, 'BtV')
    store.BtV = B' * V;
end

% Compute the products D*(V'*AV), D*(V'*MV) and D*(V'*N_iV)
if ~isfield(store, 'VtAV')
    store.VtAV = V'*store.AV;
    store.VtAV = symmetrize(store.VtAV); % Make numerically symmetric
end
if ~isfield(store, 'DVtAV')
    store.DVtAV = d .* store.VtAV;
end
if ~isfield(store, 'VtMV')
    if is_M
        store.VtMV = V'*store.MV;
        store.VtMV = symmetrize(store.VtMV); % Make numerically symmetric
    else
        store.VtMV = speye(r);
    end
end
if ~isfield(store, 'DVtMV')
    if is_M
        store.DVtMV = d .* store.VtMV;
    else
        store.DVtMV = X.D;
    end
end
if ~isfield(store, 'VtNV')
    store.VtNV = pagemtimes(V', store.NV);
    store.VtNV = symmetrize(store.VtNV); % Make numerically symmetric
end
if ~isfield(store, 'DVtNV')
    store.DVtNV = bsxfun(@times, d, store.VtNV);
end

% Compute the products C*V and V'*C*V
if ~isfield(store, 'CV')
    store.CV = B * store.BtV;
end
if ~isfield(store, 'VtCV')
    store.VtCV = V' * store.CV;
    store.VtCV = symmetrize(store.VtCV); % Make numerically symmetric
end
end
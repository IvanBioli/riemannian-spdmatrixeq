function [f, store] = lyap_psd_cost(X, A, M, N, B, store)
% lyap_psd_cost - Computes the energy-norm based cost function at
% the current iterate
%
% Syntax:
%   [f, store] = lyap_psd_cost(X, A, M, N, B, store)
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
%   - f    : Cost
%   - store: Updated store structure
%

% Precomputations
store = lyap_psd_prepare_cost(X, A, M, N, B, store);

% Compute f
f = mat_inner(store.DVtAV', store.DVtMV) - ...
    .5 * mat_inner(pagetranspose(store.DVtNV), store.DVtNV) - ...
    diag(X.D)'*diag(store.VtCV);
end


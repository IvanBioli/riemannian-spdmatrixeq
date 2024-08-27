function [g, store] = lyap_psd_grad(X, A, M, N, B, store)
% lyap_psd_grad - Computes the Riemannian gradient of the
% energy-norm based cost function.
%
% Syntax:
%   [g, store] = lyap_psd_grad(X, A, M, N, B, store)
%
% Inputs:
%   - X    : Current iterate
%   - A, M : Coefficient matrices of the multiterm linear matrix equation
%   - N    : Cell array containing the other coefficient matrices of the
%            multiterm linear matrix equation
%   - B    : Factor of the righ-hand-side B*B'
%   - store: Manopt's StoreDB struct for the current iterate
%
% Outputs:
%   - g    : Riemannian gradient
%   - store: Updated store structure
%

% Precomputations
store = lyap_psd_prepare_cost(X, A, M, N, B, store);

% Compute L(X)V
ZV = store.AV * store.DVtMV + store.MV * store.DVtAV -...
    sum(pagemtimes(store.NV, store.DVtNV), 3) - store.CV;

% Compute gradient
VtZV = X.V' * ZV; VtZV = symmetrize(VtZV);
g.M = VtZV;
g.Vp =  ZV - X.V*VtZV;

end
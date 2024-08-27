function [g, store] = lyap_psd_grad_MXM(X, A, M, N, B, store, Mass_d)
% lyap_psd_grad_MXM - Computes the Riemannian gradient
% of the energy-norm based cost function, in the metric induced by
% P(X) = Mass * X * Mass.
%
% Syntax:
%   [g, store] = lyap_psd_grad_MXM(X, A, M, N, B, store, Mass_d)
%
% Inputs:
%   - X     : Current iterate
%   - A, M  : Coefficient matrices of the multiterm linear matrix equation
%   - N     : Cell array containing the other coefficient matrices of the
%              multiterm linear matrix equation
%   - B     : Factor of the righ-hand-side B*B'
%   - Mass_d: Decomposition of the matrix Mass
%
% Outputs:
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
VtZV = X.V' * ZV;
g.M = VtZV;
g.MVp =  ZV - X.MV * VtZV;
g.Vp = Mass_d \ g.MVp;

end
function [g, store] = sylv_grad_EXD(X, A, B, C, store, E_d, D_d)
% sylv_grad_EXD - Computes the Riemannian gradient of the
% energy-norm based cost function, in the metric induced by P(X) = EXD.
%
% Syntax:
%   [g, store] = sylv_grad_EXD(X, A, B, C, store, E_d, D_d)
%
% Inputs:
%   - X    : Current iterate
%   - A, B : Cell arrays containing the coefficient matrices of the
%            multiterm linear matrix equation
%   - C    : Structure of the right-hand-side. It is a structure with
%            fields L and R, and the RHS is C.L * C.R'
%   - store: Manopt's StoreDB struct for the current iterate
%   - E_d  : Decomposition of the matrix E
%   - D_d  : Decomposition of the matrix D
%
% Outputs:
%   - g    : Riemannian gradient
%   - store: Updated store structure
%

% Precomputations
store = sylv_prepare_cost(X, A, B, C, store);
V = X.V;

% Compute Z * V and U' * Z
ZV = sum(pagemtimes(store.AU, store.SVtBtV), 3);
UtZ = sum(pagemtimes(store.UtAUS, "none", store.BV, "transpose"), 3);
ZV = ZV - C.L * store.CRtV;
UtZ = UtZ - store.CLtU' * C.R';

% Compute the factors in the gradient
g.M = UtZ * V;
g.EUp = ZV - X.EU * g.M;
g.DVp = UtZ' - X.DV * g.M';

% NOTE: g.Up and g.Vp are not needed to run optimization if using the
% approximate Newton method with B(X)=AXD+EXB. Still, the cost of
% computing them very low compared to solving the Sylvesyer equation,
% hence we compute them for compatibility with Manopt (need to compute
% gradnorm)
g.Up = E_d \ g.EUp;
g.Vp = D_d \ g.DVp;
end


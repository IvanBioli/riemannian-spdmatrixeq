function [g, store] = sylv_grad(X, A, B, C, store)
% sylv_grad - Computes the Riemannian gradient of the
% energy-norm based cost function.
%
% Syntax:
%   [g, store] = sylv_grad(X, A, B, C, store)
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
%   - g    : Riemannian gradient
%   - store: Updated store structure
%

% Precomputations
store = sylv_prepare_cost(X, A, B, C, store);
U = X.U; V = X.V;

% Compute Z * V and U' * Z
ZV = sum(pagemtimes(store.AU, store.SVtBtV), 3);
UtZ = sum(pagemtimes(store.UtAUS, "none", store.BV, "transpose"), 3);
ZV = ZV - C.L * store.CRtV;
UtZ = UtZ - store.CLtU' * C.R';

% Compute the factors in the gradient
g.M = UtZ * V;
g.Up = ZV - U * g.M;
g.Vp = UtZ' - V * g.M';
end


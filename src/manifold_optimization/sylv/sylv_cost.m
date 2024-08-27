function [f, store] = sylv_cost(X, A, B, C, store)
% sylv_cost - Computes the energy-norm based cost function at
% the current iterate
%
% Syntax:
%   [f, store] = sylv_cost(X, A, B, C, store)
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
%   - f    : Cost
%   - store: Updated store structure
%

store = sylv_prepare_cost(X, A, B, C, store);
f = 0.5 * mat_inner(store.UtAUS, store.SVtBV) - mat_inner(store.CLtU .* diag(X.S)', store.CRtV);
end


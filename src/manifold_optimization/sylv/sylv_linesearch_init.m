function [t, store] = sylv_linesearch_init(x, d, store, problem)
% sylv_linesearch_init - Computes the initial step-size for the
% linesearch procedure.
%
% Syntax:
%   [t, store] = sylv_linesearch_init(x, d, store, problem)
%
% Inputs:
%   - x      : Current iterate
%   - d      : Descent direction (tangent vector at X)
%   - store  : Manopt's StoreDB struct for the current iterate
%   - problem: Manopt's problem structure. It must have a field proj_ehess
%              which is a function computing the projected Euclidean Hessian
%
% Outputs:
%   - t    : Initial step-size
%   - store: Updated store struct
%

[Hessd, store] = problem.proj_ehess(x, d, store);
numo = problem.M.inner(x, store.grad__, d);
deno = problem.M.inner(x, d, Hessd);
t = - numo/deno;
end


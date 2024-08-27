function X = fadi_precond(X, A, B, p, q, D, E)
% fadi_precond - Applies the fADI preconditioner to X, i.e. applies fADI
% to approximately solve A*X*D' + E*X*B' = X.U*X.V'.

if isempty(p)
    return
end
if nargin < 6
    [U, V] = fadi(A, B, X.U, -X.V, p, q);
else
    [U, V] = fadi(A, B, X.U, -X.V, p, q, D, E);
end
X.U = U; X.V = V;
end


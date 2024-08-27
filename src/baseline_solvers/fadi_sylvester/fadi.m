function [W, Y] = fadi(A, B, U, V, p, q, D, E)
% fadi - Solves the Sylvester equation A*X*D' + E*X*B' + U*V' = 0.
%   Assumes that the eigenvalues of (A,E) and (B,D) have positive real part.
%   If D, E are not specified, it assumes that they are the indentity matrix,
%   i.e. solves the Sylvester equation A*X + X*B' + U*V' = 0.
%   The solution is returned in factored form X = W*Y'.
%
% Syntax:
%   [W, Y, res] = fadi(A, B, U, V, p, q, D, E)
%   [W, Y, res] = fadi(A, B, U, V, p, q)
%
% Inputs:
%   - A, B, D, E: Coefficient matrices
%   - U, V      : RHS factors
%   - p, q      : ADI poles

[m, r] = size(U);
n = size(V, 1);
steps = length(p);

W = zeros(m, steps*r);
Y = zeros(n, steps*r);

if nargin < 7
    Im = eye(m, 'like', A);
    In = eye(n, 'like', B);
    E = Im; D = In;
end

W(:, 1:r) = (A - q(1)*E) \ U;
Y(:, 1:r) = (B + p(1)*D) \ V;

k = size(U, 2);

for i = 2 : steps
    W(:, r*(i-1)+1 : r*i) = (A - q(i)*E) \ ((A - p(i-1)*E) * W(:, r*(i-2)+1 : r*(i-1)));
    Y(:, r*(i-1)+1 : r*i) = (B + p(i)*D) \ ((B + q(i-1)*D) * Y(:, r*(i-2)+1 : r*(i-1)));
end

D = kron(sqrt(diag(abs(p-q))), eye(k));
W = -W * D;
Y = Y * D;

end


function Z = fadi_lyap(A, B, p, M)
% fadi_lyap - Solves the Lyapunov equation A*X*M' + M*X*A' = B*B'.
% Assumes that the eigenvalues of (A, M) have positive real part.
% If M is not specified, it assumes that it is the indentity matrix,
% i.e. solves the Lyapunov equation A*X + X*A' = B*B'.
% The solution is returned in factored form X = Z*Z'.
%
% Syntax:
%   Z = fadi_lyap(A, B, p, M)
%
% Inputs:
%   - A, M: Coefficient matrices
%   - B   : RHS factors
%   - p   : ADI poles

[n, r] = size(B);
if nargin < 4
    In = eye(n, 'like', A);
    M = In;
end

% Preallocate
steps = length(p);
Z = zeros(n, steps*r);

% fADI iterations
Z(:, 1:r) = (A + p(1)*M) \ B;
for i = 2 : steps
    Z(:, r*(i-1)+1 : r*i) = (A + p(i)*M) \ ((A - p(i-1)*M) * Z(:, r*(i-2)+1 : r*(i-1)));
end
Z = Z * kron(sqrt(diag(2 * abs(p))), eye(r));
end


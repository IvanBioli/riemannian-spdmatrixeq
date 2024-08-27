function [X] = adi_fullrank(F, A, B, p, q, D, E)
% adi_fullrank - Solves the Sylvester equation A*X*D' + E*X*B' = F.
%   Assumes that the eigenvalues of (A,E) and (B,D) have positive real part.
%   If D, E are not specified, it assumes that they are the indentity matrix,
%   i.e. solves the Sylvester equation A*X + X*B' = F.
%
% Syntax:
%   [X] = adi_fullrank(F, A, B, p, q, D, E)
%   [X] = adi_fullrank(F, A, B, p, q)
%
% Inputs:
%   - A, B, D, E: Coefficient matrices
%   - F         : RHS
%   - p, q      : ADI poles

[m, n] = size(F);
steps = length(p);

X = zeros(m, n);

if nargin < 7
    Im = eye(m, 'like', A);
    In = eye(n, 'like', B);
    E = Im; D = In;
end

for i = 1 : steps
    % X = (A - q(i)*E) \ (F - X * (B + q(i)*D));
    % X = (F - (A - p(i)*E) * X) / (B + p(i)*D);
    X = (A - q(i)*E) \ [(p(i) - q(i)) * F + (A - p(i)*E) * X * (B + q(i)*D)] / (B + p(i) * D);
end

end


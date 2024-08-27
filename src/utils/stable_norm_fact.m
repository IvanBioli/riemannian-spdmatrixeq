function n = stable_norm_fact(X, type)
% stable_norm_fact - Computes the 2-norm/frobenius norm of a factored
% matrix, in a numerically stable way.
%   For instance, to compute the Frobenius norm of U*S*V' it first computes
%   the QR factorizations U = Q_U * R_U, V = Q_V * R_V and then exploits
%   that norm(U*S*V', 'fro') = norm(R_U*S*R_V', 'fro').
%
% Syntax:
%   [n] = stable_norm_fact(X, type)
%   [n] = stable_norm_fact(X)
%
% Inputs:
%   - X   : Struct array representing the matrix in factored form. The
%       allowed factorizations are X.U * X.V', X.U * X.S * X.V' and
%       X.V * X.D * X.V'.
%   - type: Norm type. Can be 2 (for 2-norm) or "fro" (for Frobenius norm).
%       Default is 'fro'.
%
% Outputs:
%   - n: Norm
%


if nargin < 2
    type = "fro";
else
    assert((string(type) == "fro") || (type == 2))
end
if isfield(X, 'V') && ~isfield(X, 'U')
    r = qr(full(X.V), 0);
    n = norm(r * X.D * r', type);
elseif isfield(X, 'U')
    rU = qr(full(X.U), 0);
    rV = qr(full(X.V), 0);
    if isfield(X, 'S')
        matrix = rU * X.S * rV';
    else
        matrix = rU * rV';
    end
    n = norm(matrix, type);
elseif ~isstruct(X)
    n = norm(X, type);
else
    error("Unrecognized struct for R")
end
end


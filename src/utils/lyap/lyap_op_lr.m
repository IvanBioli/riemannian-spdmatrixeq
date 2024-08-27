function Y = lyap_op_lr(A, X, M, N)
% lyap_op_lr - Returns the application of the generalized Lyapunov operator 
% to the matrix represented by the structure X. The output is also in 
% factored form.
%
% Syntax:
%   [Y] = lyap_op_lr(A, X, M, N)
%
% Inputs:
%   - A, M, N: Coefficient matrices
%   - X      : Structure representing the matrix X.V * X.D * X.V'.
%
% Outputs:
%   - Y: Result of the application of the generalized Lyapunov operator to
%       the matrix represented by X. Also Y is kept in factored format and
%       represents the matrix Y.V * Y.D * Y.V'.
%

% Get number of terms and dimensions
is_M = (nargin > 2) && ~isempty(M); 
is_N = nargin > 3;
if is_N
    n_terms = length(N);
else
    n_terms = 0;
end
r = size(X.V,2);

% Build the G factor
V = zeros(size(A,1), r * (n_terms + 2));
AV = A * X.V;
if is_M
    MV = M * X.V;
else
    MV = X.V;
end
V(:, 1:r) = (AV + MV) / sqrt(2);
V(:, r+1:2*r) = (AV - MV) / sqrt(2);
for dim = 1:n_terms
    V(:, (dim+1)*r+1:(dim+2)*r) = N{dim} * X.V;
end
Y.V = V;

% Build the D factor
Y.D = kron(spdiags([1;-1;-ones(n_terms, 1)], 0, n_terms+2, n_terms+2), X.D);

end
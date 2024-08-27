function LX = sylv_op_lr(A, B, X)
% sylv_op_lr - Returns the application of the generalized Sylvester
% operator to the matrix represented by the structure X. The output is also
% in factored form.
%
% Syntax:
%   [LX] = sylv_op_lr(A, B, X)
%
% Inputs:
%   - A, B: Coefficient matrices
%   - X   : Structure representing the matrix X.U * X.V'.
%
% Outputs:
%   - LX: Result of the application of the generalized Sylvesyer operator
%         to the matrix represented by X. Also LX is kept in factored format
%         and represents the matrix LX.U * LX.V'.
%

% Get number of terms and dimensions
n_terms = length(A);
[n, r] = size(X.U);
[m, ~] = size(X.V);

% Preallocate left and right factors
LX.U = zeros(n, r * n_terms);
LX.V = zeros(m, r * n_terms);

% Fill left and right factors
for dim = 1:n_terms
    LX.U(:, (dim-1)*r+1:dim*r) = A{dim} * X.U;
    LX.V(:, (dim-1)*r+1:dim*r) = B{dim} * X.V;
end

% If present fill the diagonal factor
if isfield(X, 'D')
    LX.D = kron(spdiags(ones(n_terms, 1), 0, n_terms, n_terms), X.D);
end

end

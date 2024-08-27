function Y = sylv_op(A, B, X)
% sylv_op - Returns the application of the generalized Sylvester operator
% to X, i.e.
%   Y = A{1} * X * B{1}' + ... + A{p} * X * B{p}'

% Get number of terms and dimensions
n_terms = length(B);

% Loop over the matrices
Y = zeros(size(X));
for dim = 1:n_terms
    Y = Y + A{dim} * X * B{dim}';
end

end



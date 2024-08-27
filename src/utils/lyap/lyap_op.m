function Y = lyap_op(A, X, M, N)
% lyap_op - Returns the application of the generalized Lyapunov operator to
% X, i.e. 
%   Y = A * X * M' + M * X * A' - N{1} * X * N{1}' - ... - N{p} * X * N{p}' 

% Get number of terms
is_M = (nargin > 2) && ~isempty(M); 
is_N = nargin > 3;
if is_N
    n_terms = length(N);
else
    n_terms = 0;
end

% Loop over the matrices
if is_M
    Y = A * X * M' + M * X * A';
else
    Y = A * X + X * A';
end
for dim = 1:n_terms
    Y = Y - N{dim} * X * N{dim}';
end

end



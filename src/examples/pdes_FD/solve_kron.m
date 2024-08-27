function [U, M, b] = solve_kron(A, B, C, data)
% solve_kron - Solves the Kronecker formulation of the matrix equation from
% the finite difference discretization of a 2D PDEs with separable 
% variables and Dirichlet boundary condition
%
% Syntax:
%   [U, M, b] = solve_kron(A, B, C, data)
%
% Description:
%   U is the solution to vectorization of the linear matrix equation:
%       A{1}*U*B{1}' + ... + A{p}*U*B{p}' = C.L * C.R'  
%
% Inputs:
%   - A, B: Cell arrays containing the coefficient matrices
%   - C   : Structure representing the RHS F = C.L * C.R'
%   - data: Structure containing the Dirichlet boundary data. The field
%           data.g must be a two-variable function (the boundary data). If
%           there is no field data.g, it is assumed to be zero.
%
% Outputs:
%   - U: Approximate solution
%   - M: Matrix of the vectorized system
%   - b: RHS of the vectorized system
%

% Get problem dimension and allocate space for M and b
n = size(A{1}, 1);
M = sparse(n^2, n^2);
b = sparse(n^2, 1);

% Build RHS of the system
for i = 1:length(A)
    M = M + kron(B{i}, A{i});
end

% Build LHS of the system
for i = 1:size(C.L,2)
    b = b + kron(C.R(:, i), C.L(:, i));
end

% Solve and reshape
u_vec = M\b; 
u_mat = reshape(u_vec, [n, n]); 
U = zeros(n+2,n+2); 
U(2:end-1, 2:end-1) = u_mat; 

% Impose Dirichlet boundary conditions
if isfield(data, 'g')
    x = linspace(0, 1, n+2);
    U(1, :) = data.g(0, x);
    U(end, :) = data.g(1, x);
    U(:, 1) = data.g(x, 0);
    U(:, end) = data.g(x, 1);
end

end


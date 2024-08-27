function [A, B, C] = assemble_pdes_FD(n, data)
% assemble_pdes_FD - Assembles the multiterm linear matrix equation
% obtained from the finite difference discretization of a 2D elliptic PDE
% with separable variables.
%   The PDE is of the form
%       -div(k(x,y) * grad(u(x,y)) = f(x,y) in ]0,1[^2                (*)
%                           u(x,y) = g(x,y) in boundary([0,1]^2)      (*)
%   where
%       k(x,y) = eps + k_coeff(1) * k_x{1}(x) * k_y{1}(y) + ... +
%                           k_coeff(p) * k_x{p}(x) *k_y{p}(y)         (**)
%       f(x,y) = f_x{1}(x) * f_y{1}(y) + ... + f_x{p}(x) * f_y{p}(y)  (**)
%   For the resulting diffusion multiterm linear matrix equation we use the
%   notation
%           A{1} * X * B{1}' + ... + A{ell} * X * B{ell}' = F
%
% Syntax:
%   [A, B, C] = assemble_pdes_FD(n, data)
%
% Inputs:
%   - n   : Number of mesh-points per dimension. The mesh width is h = 1/n.
%   - data: Structure describing the diffusion, source and boundary data of
%       equation (*). Using the same notation as above, data is a
%       structure with fields data.dataname where dataname is is one of
%       the following:
%           - k_x, k_y: Cell arrays of functions defining the diffusion term
%           - k: function of diffusion term k(x,y) in (**)
%           - eps: Scalar
%           - k_coeff: Array
%           - f_x, f_y: Cell arrays of functions defining the source term
%           - g: Function defining the boundary data
%
% Outputs:
%   - A, B: Cell arrays containing the coefficient matrices of the
%           multiterm linear matrix equation
%   - C   : Structure representing the RHS F = C.L * C.R'
%

% Merge with default data
data = merge_default_data(data);

% Precomputations
h = 1/n;
x = linspace(0, 1, n+1)';
x_half = 0.5 * (x(1:end-1) + x(2:end));
x_minus = x_half(2:end); x_plus = x_half(1:end-1);
x_int = x(2:end-1);
e = ones(n-1,1);

% Assemble LHS of the system
L = spdiags([-e 2*e -e], -1:1, n-1, n-1);
I = speye(n-1);
A{1} = I; A{2} = data.eps * L;
B{1} = data.eps * L; B{2} = I;
for i = 1:length(data.k_x)
    % Define discretization matrices
    A_x = spdiags(data.k_x{i}(x_int), 0, n-1, n-1);
    B_x = spdiags([...
        -data.k_x{i}(x_minus),...
        data.k_x{i}(x_plus)+data.k_x{i}(x_minus),...
        -data.k_x{i}(x_plus)...
        ], -1:1, n-1, n-1);
    A_y = spdiags(data.k_y{i}(x_int), 0, n-1, n-1);
    B_y = spdiags([...
        -data.k_y{i}(x_minus),...
        data.k_y{i}(x_plus)+data.k_y{i}(x_minus),...
        -data.k_y{i}(x_plus)...
        ], -1:1, n-1, n-1);
    
    % Assign to cell arrays
    c = sqrt(data.k_coeff(i));
    A{2*i+1} = c * B_x;
    A{2*i+2} = c * A_x;
    B{2*i+1} = c * A_y;
    B{2*i+2} = c * B_y;
end

% Assemble RHS of the system
C = struct(); C.L = []; C.R = [];
if isfield(data, 'f_x')
    for i = 1:length(data.f_x)
        C.L = [C.L, data.f_x{i}(x_int)];
        C.R = [C.R, data.f_y{i}(x_int)];
    end
    C.L = h * C.L; C.R = h * C.R;
end
% Add contribute from Dirichlet BC
if isfield(data, 'g')
    vl = sparse(1, 1, 1, n-1, 1);
    vr = sparse(n-1, 1, 1, n-1, 1);
    % Left boundary
    vec = data.k(h/2, x_int) .* data.g(0, x_int);
    C.L = [C.L, vl]; C.R = [C.R, vec];
    % Right boundary
    vec = data.k(1-h/2, x_int) .* data.g(1, x_int);
    C.L = [C.L, vr]; C.R = [C.R, vec];
    % Upper boundary
    vec = data.k(x_int, 1-h/2) .* data.g(x_int, 1);
    C.L = [C.L, vec]; C.R = [C.R, vr];
    % Lower boundary
    vec = data.k(x_int, h/2) .* data.g(x_int, 0);
    C.L = [C.L, vec]; C.R = [C.R, vl];
end

end


function data_new = merge_default_data(data)
% merge_default_data - Merges the data structure given by the user with
%   the default option structure. User-defined options have precedence
%   over the default ones.
%
% Syntax:
%   [data_new] = merge_default_data(data)
%
% Inputs:
%   - options: Options structure with field options.fieldname, where
%       field is one of the following and the default value is indicated
%       between parentheses:
%           - eps (1): see equations (*) and (**)
%           - k_coeff (ones(length(data.k_x), 1);): see equations (*) and
%             (**)
%
% Outputs:
%   - data_new: Data structure containing all the fields above, with
%     default values in fields that have not been provided by the user
%

% Define default data
data_default.eps = 1;
data_default.k_coeff = ones(length(data.k_x), 1);

% Merge data with defaults
data_new = mergeOptions(data_default, data);
end
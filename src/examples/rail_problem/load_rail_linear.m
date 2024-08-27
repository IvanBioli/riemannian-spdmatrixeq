% Description: Reachability Gramian of a linear control system (not included
% in the report)
%   Loads the matrices and selects the first column of the RHS matrix B

% Load matrices
if ~exist("k_rail", "var")
    k_rail = 0;
end
eqn = mess_get_linear_rail(k_rail);
n = size(eqn.A_, 1);
A = -eqn.A_; % A is SPD in our formulation
M = eqn.E_;
B = eqn.B;
N = {}; % The last matrix is zero

% Rescale matrices
scaling = 1;
A = A * scaling;
M = M * scaling;
B = B * scaling;

% Compute the reordering and Cholesky factorizations
start = tic();
p = symamd(A);
A_old = A; A = A(p,p);
M_old = M; M = M(p, p);
B = B(p, 1);
M_d = decomposition(M);
time_precomp = toc(start);

% % Plot sparsity pattern
% set(groot,'defaulttextinterpreter','latex');
% set(groot, 'defaultAxesTickLabelInterpreter','latex');
% set(groot, 'defaultLegendInterpreter','latex');
% fig = figure();
% t = tiledlayout(2, 3,'TileSpacing','Compact');
% nexttile; spy(A_old);
% nexttile; spy(A);
% nexttile; spy(chol(A));
% nexttile; spy(M_old);
% nexttile; spy(M);
% nexttile; spy(chol(M));
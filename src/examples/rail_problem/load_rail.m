% Description: Reachability Gramian of a bilinear control system
%   Loads the matrices

% Load matrices
if ~exist("k_rail", "var")
    k_rail = 0;
end
eqn = mess_get_bilinear_rail_modified(k_rail);
n = size(eqn.A_, 1);
A = -eqn.A_; % A is SPD in our formulation
M = eqn.E_;
B = eqn.B;
N = eqn.N_(1:6); % The last matrix is zero

% Rescale matrices
scaling = 1e3;
A = A * scaling;
M = M * scaling;
B = B * scaling;
for i = 1:6
    N{i} = N{i} * scaling;
end

% Compute the reordering and Cholesky factorizations
start = tic();
p = symamd(A);
A_old = A; A = A(p,p);
M_old = M; M = M(p, p);
M_d = decomposition(M);
for i = 1:length(N)
    N_curr = N{i};
    N{i} = N_curr(p,p);
end
B = B(p,:);
B = B(:, [1,7]);
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
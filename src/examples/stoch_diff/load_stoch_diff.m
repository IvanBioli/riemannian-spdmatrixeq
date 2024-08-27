% Description: Stochastic Galerkin Matrix Equations
%   Loads the matrices

% Load datafiles
if TP == 5
    load examples/stoch_diff/MultiRB/TP_five.mat;
    load ("examples/stoch_diff/KL_data_TP5.mat", "aa");
    if ~exist("ranks", "var")
    end
elseif TP == 2
    load examples/stoch_diff/MultiRB/TP_two.mat;
    load ("examples/stoch_diff/KL_data_TP2.mat", "aa");
    if ~exist("ranks", "var")
    end
else
    error("Unknown test problem")
end

% Rewrite the equation compatibly with our formulation
m = size(K{1}, 1); n = size(G{1}, 1);
F = reshape(fnew, [m, n]); C.L = F(:, 1); C.R = sparse(1, 1, 1, n, 1);
assert(norm(F - C.L*C.R', 'fro') < eps)
ell = length(K);
A = cell(ell, 1); B = cell(ell, 1);
for i = 1:ell
    A{i} = 0.5*(K{i}+K{i}');
    B{i} = G{i};
    assert(norm(K{i} - A{i}, 'fro') / norm(K{i}, 'fro') < eps)
    assert(issymmetric(G{i}));
end

% Define the inverse of Amean as an operator
Amean_old = Amean;
Amean = 0.5 * (Amean + Amean');
assert(norm(Amean - Amean_old, 'fro') / norm(Amean_old, 'fro') < eps)
% Permute the system based on Amean
start = tic();
pa = dissect(Amean);
Amean = Amean(pa, pa);
Amean_d = decomposition(Amean);
time_precomp_A = toc(start);
for i = 1:length(A)
    A{i} = A{i}(pa,pa);
end
C.L = C.L(pa,:);

% Compute the coefficients for Gmean
% deno = norm(A{1}, 'fro')^2;
% a_bar = zeros(length(B), 1);
% a_bar(1) = 1;
% for i = 2:length(B)
%     a_bar(i) = mat_inner(A{i}, A{1}) / deno;
% end
a_bar = mean(aa, 1);

% Permute the system based on Gmean
Gmean = sparse(size(G{1}, 1), size(G{1}, 2));
for i = 1:ell
    Gmean = Gmean + a_bar(i) * B{i};
end
Gmean_old = Gmean;
start = tic();
pb = dissect(Gmean);
Gmean = Gmean(pb, pb);
Gmean_d = decomposition(Gmean);
time_precomp_G = toc(start);
for i = 1:length(B)
    B_curr = B{i};
    B{i} = B_curr(pb,pb);
end
C.R = C.R(pb,:);

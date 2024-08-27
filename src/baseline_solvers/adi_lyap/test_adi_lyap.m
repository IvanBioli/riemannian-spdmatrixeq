% Description: Test script for fadi_lyap and adi_lyap.m

clear all; clc;
%% Test of fadi_lyap
fprintf("fadi_lyap:\n")
n = 100;
npoles = 100;
A = randn(n); A = A * A';

eA = eig(A); a = min(eA); b = max(eA);
p = zolotarev_poles(npoles, a, b);
k = 2;
B = randn(n, k);
Z = fadi_lyap(A, B, p);
X = Z * Z';
res = A * X + X * A' - B * B';
fprintf("\t without M: %e\n", norm(res, 'fro') / norm(B * B', 'fro'))

M = randn(n); M = M * M';
eAM = eig(A, M); a = min(eAM); b = max(eAM);
p = zolotarev_poles(npoles, a, b);
k = 2;
B = randn(n, k);
Z = fadi_lyap(A, B, p, M);
X = Z * Z';
res = A * X * M' + M' * X * A' - B * B';
fprintf("\t with M: %e\n", norm(res, 'fro') / norm(B * B', 'fro'))

%% Test of adi_lyap
clear all;
fprintf("adi_lyap:\n")
[A, N, B] = assemble_heat_conduction(10, [1], [], 1); n = 100;
npoles = 1000;
[Q, ~] = qr(rand(100)); M = Q * diag(1 + rand(100, 1)) * Q';
eAM = eig(full(A), M); a = min(eAM); b = max(eAM);
p = zolotarev_poles(npoles, a, b);

k = 2;
X.V = randn(n, k); X.D = diag(randn(k, 1));
Y = adi_lyap(A, X, p, struct(), M, {});
X_ = X.V * X.D * X.V';
Y_ = Y.V * Y.D * Y.V';
res = A * Y_ * M' + M * Y_ * A' - X_;
fprintf("\t without N: %e\n", norm(res, 'fro') / norm(X_, 'fro'))

options.trunc_tol = 1e-16;
Y = adi_lyap(A, X, p, options, M, N);
Y_ = Y.V * Y.D * Y.V';
res = A * Y_ * M' + M * Y_ * A' - X_;
for i = 1:length(N)
    res = res - N{i} * Y_ * N{i}';
end
fprintf("\t with N: %e\n", norm(res, 'fro') / norm(X_, 'fro'))

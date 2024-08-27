% Description: Test script for fadi.m

n = 1024;
npoles = 100;

A = randn(n);
A = A * A';
B = randn(n);
B = B * B';
D = randn(n);
D = D * D';
E = randn(n);
E = E * E';

eigAE = eig(A, E); a = min(eigAE); b = max(eigAE);
eigBD = eig(B, D); c = min(eigBD); d = max(eigBD);

[p,q] = zolotarev_poles(npoles, a, b, c, d);

U = randn(n, 2);
V = randn(n, 2);

[W, Y] = fadi(A, B, U, V, p, q, D, E);

X = W * Y';
norm(A*X*D + E*X*B' + U*V', 'fro') / norm(U*V', 'fro')
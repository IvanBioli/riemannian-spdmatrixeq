% Description: Test script for the example "Finite difference discretization of 2D PDEs
% with separable variables"

clear all;
% %% PROBLEM DEFINITION
% syms x y a
% k(x,y) = 1 + a * cos(2 * pi * x);
% u(x, y) = sin(2*pi*x) * y * (1-y);
% f = diff(k * diff(u, 1, x), 1, x) + diff(k * diff(u, 1, y), 1, y);
% simplify(f)
%
% a = 0.5;
% data.k =@(x,y) 1 + a * cos(2 * pi* x) .* ones(size(y));
% data.k_x{1} =@(x) a * cos(2 * pi *x);
% data.k_y{1} =@(y) ones(size(y));
%
% data.f_x{1} =@(x) sin(2*pi*x).*(a*cos(2*pi*x) + 1);
% data.f_y{1} =@(y) 4*y*pi^2 .* (y-1);
% data.f_x{2} =@(x) -2*sin(2*pi*x).*(a*cos(2*pi*x) + 1);
% data.f_y{2} =@(y) ones(size(y));
% data.f_x{3} =@(x) cos(2*pi*x).*sin(2*pi*x);
% data.f_y{3} =@(y) 4*a*y*pi^2.*(y - 1);
% u =@(x,y) sin(2*pi*x) .* y .* (1-y);
% data.g = u;

syms x y a
k(x,y) = 1 + a * sin(x) + cos(y);
u(x, y) = cos(x) + sin(y);
f = -diff(k * diff(u, 1, x), 1, x) - diff(k * diff(u, 1, y), 1, y);
simplify(f)

%% PROBLEM DEFINITION
a = 1;
data.k =@(x,y) 1 + a * sin(x) * ones(size(y)) + cos(y) * ones(size(x)) ;
data.k_x{1} =@(x) a * sin(x);
data.k_y{1} =@(y) ones(size(y));
data.k_x{2} =@(x) ones(size(x));
data.k_y{2} =@(y) cos(y);

data.f_x{1} =@(x) cos(x);
data.f_y{1} =@(y) 1 + cos(y);
data.f_x{2} =@(x) 1 + a * sin(x);
data.f_y{2} =@(y) sin(y);
data.f_x{3} =@(x) 2 * a * cos(x) .* sin(x);
data.f_y{3} =@(y) ones(size(y));
data.f_x{4} =@(x) ones(size(x));
data.f_y{4} =@(y) 2 * cos(y) .* sin(y);

u =@(x,y) cos(x) + sin(y);
data.g = u;

%% SOLVE AND PLOT THE ERROR
n_vec = 2.^(2:8);
err_vec = [];
for n = n_vec
    [Y, X] = meshgrid(linspace(0, 1, n+1)); %To me U_ij = u(ih, jh) --> flip Y and X
    [A, B, C] = assemble_pdes_FD(n, data);
    U = solve_kron(A, B, C, data);
    U_exact = u(X, Y);
    err = max(abs(U-U_exact),[], 'all') / max(abs(U_exact),[], 'all');
    err_vec = [err_vec, err];
end

figure()
loglog(n_vec, err_vec, '.-', 'DisplayName', 'Error')
hold on
loglog(n_vec, 1./n_vec.^2, '--', 'DisplayName', '1/n^2')
legend()
% Description: Test script for the example "Finite difference discretization of 2D PDEs
% with separable variables"

clearvars -except n;
%% ASSEMBLE PROBLEM
n_terms = 3;

% Gather informations
% n = 1000;
a_vec = [1, 10];
cond_vec = zeros(size(a_vec)); sigma_vec = {};
for i = 1:length(a_vec)
    a = a_vec(i);
    [Y, X] = meshgrid(linspace(0, 1, n+1)); %To me U_ij = u(ih, jh) --> flip Y and X
    [data, u] = example1(a, n_terms);
    [A, B, C] = assemble_pdes_FD(n, data);
    [U, M, b] = solve_kron(A, B, C, data);
    figure();
    surf(X, Y, U);
    cond_vec(i) = condest(M);
    sigma_vec{i} = svd(U);
end
%%
figure()
loglog(a_vec, cond_vec, '.-', 'DisplayName', 'cond(kron(L)))')
legend()
title("Conditioning number")
xlabel('a')

fig = figure();
for i = 1:length(a_vec)
    s = sigma_vec{i};
    y = sqrt(cumsum(s.^2, 'reverse') / sum(s.^2, 'all'));
    semilogy(y(y>1e-15), 'DisplayName', "a = "+num2str(a_vec(i)))
    hold on
end
xlabel('i')
ylabel('$\sqrt{\sum_{j=i+1}^n\sigma_j^2 / \sum_{j=1}^n\sigma_j^2}$', 'Interpreter','latex')
legend()
title("Relative error by low-rank approximation (n = " + num2str(n) +")")
exportgraphics(fig, "../report/figures/pdes_FD_lowrankapprox_n" + num2str(n) + ".pdf")
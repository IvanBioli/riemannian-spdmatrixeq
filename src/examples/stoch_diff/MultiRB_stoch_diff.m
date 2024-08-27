% Description: Stochastic Galerkin Matrix Equations
%   Executes experiments with MultiRB solver

clearvars -except seed TP tol p_option_vec;
if ~exist("p_option_vec", "var")
    p_option_vec = [1, 2];
end
results = {};
names = {};
for p_option = p_option_vec
    SGFEM_matdriver;
    results{p_option} = info;
    if p_option == 1
        names{p_option} = "MultiRB parameter-free";
    else
        names{p_option} = "MultiRB parameter-dependent";
    end
    clearvars -except seed p_option TP tol results names
end
save("results/MultiRB_TP" + num2str(TP) + ".mat", "results", "names", "tol")
% Description: Stochastic Galerkin Matrix Equations
%   Executes experiments with Riemannian optimization

clearvars -except seed TP fixed_rank ...
    sanity_checks whatruns_firstord whatruns_rtr ...
    firstord_only rtr_only rank_adaptivity opts_adaptivity...
    plot_and_save rel_tolgradnorm; clc;
%% GENERATING PROBLEM INSTANCE
% Load data
% TP = 5;
load_stoch_diff;

% General options
rhs_norm = stable_norm_fact(struct('U', C.L, 'V', C.R));
if exist("rel_tolgradnorm", "var")
    tolgradnorm = rel_tolgradnorm * rhs_norm;
else
    if TP == 5
        tolres = 1e-6 * rhs_norm;
    elseif TP == 2
        tolres = 1e-5 * rhs_norm;
    end
    tolgradnorm = tolres / 2;
end
statsfun_ = @(problem, x, stats, store) sylv_statsfun(problem, x, stats, store);

% Options for conjugategradient
% RGD
opts_rgd = struct();
opts_rgd.beta_type = "steep";
opts_rgd.tolgradnorm = tolgradnorm;
opts_rgd.minstepsize = 1e-4 * opts_rgd.tolgradnorm;
opts_rgd.maxiter = 100;
opts_rgd.statsfun = statsfun_;
% RCG
opts_rcg = struct();
opts_rcg.beta_type = "H-S-SATO";
opts_rcg.tolgradnorm = tolgradnorm;
opts_rcg.minstepsize = 1e-4 * opts_rcg.tolgradnorm;
opts_rcg.maxiter = opts_rgd.maxiter;
opts_rcg.statsfun = statsfun_;
cell_optsfirstord = {opts_rgd, opts_rcg};

% Options for trustregions
opts_rtr_rhess = struct();
opts_rtr_rhess.projehess = false;
opts_rtr_rhess.tolgradnorm = tolgradnorm;
opts_rtr_rhess.maxiter = 100;
opts_rtr_rhess.maxinner = 250;
opts_rtr_rhess.statsfun = statsfun_;
opts_rtr_projehess = opts_rtr_rhess;
opts_rtr_projehess.projehess = true;
cell_optsrtr = {opts_rtr_rhess, opts_rtr_projehess};

% Set experiments to run
if ~exist("rank_adaptivity", "var")
    rank_adaptivity = false;
end
if rank_adaptivity
    name_prefix = "RRAM + ";
    name_postfix = "";
else
    name_prefix = "";
    name_postfix = ", $r=" + num2str(fixed_rank) + "$";
end

if rank_adaptivity
    fixed_rank = opts_adaptivity.r_start;
end

if ~exist("whatruns_firstord", "var")
    whatruns_firstord = ones(1, 7);
end
if ~exist("whatruns_rtr", "var")
    whatruns_rtr = ones(1, 5);
end
whatruns_final = [];

%% OPTIMIZATION

% DEFINE THE MANIFOLD, THE COST AND THE GRADIENT
% Problem with standard metric
problem.A = A; problem.B = B; problem.C = C;
problem.M = fixedrankembeddedfactory(m, n, fixed_rank);
problem.cost = @(X, store) sylv_cost(X, A, B, C, store);
problem.grad = @(X, store) sylv_grad(X, A, B, C, store);
problem.proj_ehess = @(X, eta, store) sylv_hess(X, eta, A, B, C, store, struct('projehess', true));
problem.rhs_norm = rhs_norm;
problem.linesearch = @(x, d, store) sylv_linesearch_init(x, d, store, problem);
problem_pcg = problem;
M_handle =@(r) fixedrankembeddedfactory(m, n, r);

% Problem with metric induced by P(X) = K_mean * X
I = speye(size(B{1}, 1));
problem_changed_metric1.A = A; problem_changed_metric1.B = B; problem_changed_metric1.C = C;
problem_changed_metric1.M = fixedrank_metricEXD_factory(Amean, I, fixed_rank);
problem_changed_metric1.cost = @(X, store) sylv_cost(X, A, B, C, store);
problem_changed_metric1.grad = @(X, store) sylv_grad_EXD(X, A, B, C, store, Amean_d, I);
problem_changed_metric1.rhs_norm = rhs_norm;
problem_changed_metric1.proj_ehess = @(X, eta, store) sylv_hess_EXD(X, eta, A, B, C, store, Amean_d, I, struct('projehess', true));
problem_changed_metric1.linesearch = @(x, d, store) sylv_linesearch_init(x, d, store, problem_changed_metric1);
Mchangedmetric1_handle =@(r) fixedrank_metricEXD_factory(Amean, I, r);

% Problem with metric induced by P(X) = K_mean * X * Gmean
problem_changed_metric2.A = A; problem_changed_metric2.B = B; problem_changed_metric2.C = C;
problem_changed_metric2.M = fixedrank_metricEXD_factory(Amean, Gmean, fixed_rank, chol(Amean), chol(Gmean));
problem_changed_metric2.cost = @(X, store) sylv_cost(X, A, B, C, store);
problem_changed_metric2.grad = @(X, store) sylv_grad_EXD(X, A, B, C, store, Amean_d, Gmean_d);
problem_changed_metric2.rhs_norm = rhs_norm;
problem_changed_metric2.proj_ehess = @(X, eta, store) sylv_hess_EXD(X, eta, A, B, C, store, Amean_d, Gmean_d, struct('projehess', true));
problem_changed_metric2.linesearch = @(x, d, store) sylv_linesearch_init(x, d, store, problem_changed_metric2);
Mchangedmetric2_handle =@(r) fixedrank_metricEXD_factory(Amean, Gmean, r);

% Initialization point
init_norm = 1;
x0 = problem.M.rand(); x0.S = init_norm * x0.S;
x0_changed_metric1 = problem_changed_metric1.M.SVD2repr(x0);
x0_changed_metric2 = problem_changed_metric2.M.SVD2repr(x0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RCG/RCG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:length(cell_optsfirstord)
    opts_firstord = cell_optsfirstord{j};
    idx_shift = (j - 1) * length(whatruns_firstord);
    % Get if RCG or RGD
    if opts_firstord.beta_type == "steep"|| opts_firstord.beta_type == "S-D"
        cg_type = "R-GD";
    else
        cg_type = "R-NLCG";
    end
    if exist("firstord_only", "var") && j ~= firstord_only
        whatruns_firstord_ = zeros(size(whatruns_firstord));
    else
        whatruns_firstord_ = whatruns_firstord;
    end
    whatruns_final = [whatruns_final, whatruns_firstord_];
    
    % R-GD/R-NLCG
    if whatruns_firstord_(1)
        if ~rank_adaptivity
            [x, ~, info, ~] = conjugategradient(problem, x0, opts_firstord);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = opts_firstord;
            [~, info] = sylv_rankadaptivity(problem, @conjugategradient, M_handle, opts_adapt);
        end
        results{1+idx_shift} = info;
        names{1+idx_shift} = name_prefix + cg_type + name_postfix;
    end
    
    % Riemannian Truncated Preconditioned Richardson P(X) = K_mean * X
    if whatruns_firstord_(2)
        problem_pcg.precon = @(X, H, store) sylv_richardson_EX(X, A, B, C, Amean_d, store);
        if ~rank_adaptivity
            [x, ~, info, ~] = conjugategradient(problem_pcg, x0, opts_firstord);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = opts_firstord;
            [~, info] = sylv_rankadaptivity(problem, @conjugategradient, M_handle, opts_adapt);
        end
        info = add_time(info, time_precomp_A);
        results{2+idx_shift} = info;
        names{2+idx_shift} = name_prefix + cg_type + " + Prec. Rich. $\mathcal{P}^{(1)}X = K_0X$" + name_postfix;
    end
    
    % Riemannian Truncated Preconditioned Richardson P(X) = K * X * G_means
    if whatruns_firstord_(3)
        problem_pcg.precon = @(X, H, store) sylv_richardson_EXD(X, A, B, C, Amean_d, Gmean_d, store);
        if ~rank_adaptivity
            [x, ~, info, ~] = conjugategradient(problem_pcg, x0, opts_firstord);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = opts_firstord;
            [~, info] = sylv_rankadaptivity(problem_pcg, @conjugategradient, M_handle, opts_adapt);
        end
        info = add_time(info, time_precomp_A + time_precomp_G);
        results{3+idx_shift} = info;
        names{3+idx_shift} = name_prefix + cg_type + " + Prec. Rich. $\mathcal{P}^{(2)}X = K_0XG$" + name_postfix;
    end
    
    % R-GD/R-NLCG + Riemannian preconditioner P(X) = K_mean * X
    if whatruns_firstord_(4)
        problem_pcg.precon = @(X, H, store) sylv_precon_EX(X, H, Amean, Amean_d, store);
        if ~rank_adaptivity
            [x, ~, info, ~] = conjugategradient(problem_pcg, x0, opts_firstord);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = opts_firstord;
            [~, info] = sylv_rankadaptivity(problem_pcg, @conjugategradient, M_handle, opts_adapt);
        end
        info = add_time(info, time_precomp_A);
        results{4+idx_shift} = info;
        names{4+idx_shift} = name_prefix + cg_type + " + $\mathcal{P}^{(1)}X = K_0X$" + name_postfix;
    end
    
    % R-GD/R-NLCG + Riemannian preconditioner P(X) = K * X * G_mean
    if whatruns_firstord_(5)
        problem_pcg.precon = @(X, H, store) sylv_precon_EXD(X, H, Amean, Amean_d, Gmean, Gmean_d, store);
        if ~rank_adaptivity
            [x, ~, info, ~] = conjugategradient(problem_pcg, x0, opts_firstord);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = opts_firstord;
            [~, info] = sylv_rankadaptivity(problem_pcg, @conjugategradient, M_handle, opts_adapt);
        end
        info = add_time(info, time_precomp_A + time_precomp_G);
        results{5+idx_shift} = info;
        names{5+idx_shift} = name_prefix + cg_type + " + $\mathcal{P}^{(2)}X = K_0XG$" + name_postfix;
    end
    
    % R-GD/R-NLCG + Tools metric P(X) = K_mean * X
    if whatruns_firstord_(6)
        opts_solver = opts_firstord;
        opts_solver.tolgradnorm = opts_solver.tolgradnorm / sqrt(normest(problem_changed_metric1.M.E, 1e-4) * normest(problem_changed_metric1.M.D, 1e-4));
        opts_solver.minstepsize = 1e-4 * opts_solver.tolgradnorm;
        if ~rank_adaptivity
            [x, ~, info, ~] = conjugategradient(problem_changed_metric1, x0_changed_metric1, opts_solver);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0_changed_metric1;
            opts_adapt.opts_solver = opts_solver;
            [~, info] = sylv_rankadaptivity(problem_changed_metric1, @conjugategradient, Mchangedmetric1_handle, opts_adapt, I, Amean_d);
        end
        info = add_time(info, time_precomp_A);
        results{6+idx_shift} = info;
        names{6+idx_shift} = name_prefix + cg_type + " + Tools metric $\mathcal{P}^{(1)}X = K_0X$" + name_postfix;
    end
    
    % R-GD/R-NLCG + Tools metric P(X) = K * X * G_mean
    if whatruns_firstord_(7)
        opts_solver = opts_firstord;
        opts_solver.tolgradnorm = opts_solver.tolgradnorm / sqrt(normest(problem_changed_metric2.M.E, 1e-4) * normest(problem_changed_metric2.M.D, 1e-4));
        opts_solver.minstepsize = 1e-4 * opts_solver.tolgradnorm;
        if ~rank_adaptivity
            [x, ~, info, ~] = conjugategradient(problem_changed_metric2, x0_changed_metric2, opts_solver);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0_changed_metric2;
            opts_adapt.opts_solver = opts_solver;
            [~, info] = sylv_rankadaptivity(problem_changed_metric2, @conjugategradient, Mchangedmetric2_handle, opts_adapt, Gmean_d, Amean_d);
        end
        info = add_time(info, time_precomp_A + time_precomp_G);
        results{7+idx_shift} = info;
        names{7+idx_shift} = name_prefix + cg_type + " + Tools metric $\mathcal{P}^{(2)}X = K_0XG$" + name_postfix;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RTR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:length(cell_optsrtr)
    options_rtr = cell_optsrtr{j};
    idx_shift = length(cell_optsfirstord) * length(whatruns_firstord) + ...
        (j - 1) * length(whatruns_rtr);
    if exist("rtr_only", "var") && j ~= rtr_only
        whatruns_rtr_ = zeros(size(whatruns_rtr));
    else
        whatruns_rtr_ = whatruns_rtr;
    end
    whatruns_final = [whatruns_final, whatruns_rtr_];
    
    % Get if Riemannian Hessian or projected Euclidean Hessian
    if options_rtr.projehess
        rtr_type = "R-TR Proj. EHess";
    else
        rtr_type = "R-TR Riemann. Hess";
    end
    options_hess = struct(); options_hess.projehess = options_rtr.projehess;
    problem.hess = @(X, eta, store) sylv_hess(X, eta, A, B, C, store, options_hess);
    problem_pcg = problem;
    problem_changed_metric1.hess = @(X, eta, store) sylv_hess_EXD(X, eta, A, B, C, store, Amean_d, I, options_hess);
    problem_changed_metric2.hess = @(X, eta, store) sylv_hess_EXD(X, eta, A, B, C, store, Amean_d, Gmean_d, options_hess);
    
    % RTR
    if whatruns_rtr_(1)
        if ~rank_adaptivity
            [x, ~, info, ~] = trustregions(problem, x0, options_rtr);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = options_rtr;
            [~, info] = sylv_rankadaptivity(problem, @trustregions, M_handle, opts_adapt);
        end
        results{1+idx_shift} = info;
        names{1+idx_shift} = name_prefix + rtr_type + name_postfix;
    end
    
    % RTR + tCG + preconditioner P(X) = K_mean * X
    if whatruns_rtr_(2)
        problem_pcg.precon = @(X, H, store) sylv_precon_EX(X, H, Amean, Amean_d, store);
        if ~rank_adaptivity
            [x, ~, info, ~] = trustregions(problem_pcg, x0, options_rtr);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = options_rtr;
            [~, info] = sylv_rankadaptivity(problem_pcg, @trustregions, M_handle, opts_adapt);
        end
        info = add_time(info, time_precomp_A);
        results{2+idx_shift} = info;
        names{2+idx_shift} = name_prefix + rtr_type + " + $\mathcal{P}^{(1)}X = K_0X$" + name_postfix;
    end
    
    % RTR + tCG + preconditioner P(X) = K * X * G_mean
    if whatruns_rtr_(3)
        problem_pcg.precon = @(X, H, store) sylv_precon_EXD(X, H, Amean, Amean_d, Gmean, Gmean_d, store);
        if ~rank_adaptivity
            [x, ~, info, ~] = trustregions(problem_pcg, x0, options_rtr);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = options_rtr;
            [~, info] = sylv_rankadaptivity(problem_pcg, @trustregions, M_handle, opts_adapt);
        end
        info = add_time(info, time_precomp_A + time_precomp_G);
        results{3+idx_shift} = info;
        names{3+idx_shift} = name_prefix + rtr_type + " + $\mathcal{P}^{(2)}X = K_0XG$" + name_postfix;
    end
    
    % RTR + Tools metric P(X) = K_mean * X
    if whatruns_rtr_(4)
        opts_solver = options_rtr;
        opts_solver.tolgradnorm = opts_solver.tolgradnorm / sqrt(normest(problem_changed_metric1.M.E, 1e-4) * normest(problem_changed_metric1.M.D, 1e-4));
        opts_solver.minstepsize = 1e-4 * opts_solver.tolgradnorm;
        if ~rank_adaptivity
            [x, ~, info, ~] = trustregions(problem_changed_metric1, x0_changed_metric1, opts_solver);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0_changed_metric1;
            opts_adapt.opts_solver = opts_solver;
            [~, info] = sylv_rankadaptivity(problem_changed_metric1, @trustregions, Mchangedmetric1_handle, opts_adapt, I, Amean_d);
        end
        info = add_time(info, time_precomp_A);
        results{4+idx_shift} = info;
        names{4+idx_shift} = name_prefix + rtr_type + " + Tools metric $\mathcal{P}^{(1)}X = K_0X$" + name_postfix;
    end
    
    % RTR + Tools metric P(X) = K * X * G_mean
    if whatruns_rtr_(5)
        opts_solver = options_rtr;
        opts_solver.tolgradnorm = opts_solver.tolgradnorm / sqrt(normest(problem_changed_metric2.M.E, 1e-4) * normest(problem_changed_metric2.M.D, 1e-4));
        opts_solver.minstepsize = 1e-4 * opts_solver.tolgradnorm;
        if ~rank_adaptivity
            [x, ~, info, ~] = trustregions(problem_changed_metric2, x0_changed_metric2, opts_solver);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0_changed_metric2;
            opts_adapt.opts_solver = opts_solver;
            [~, info] = sylv_rankadaptivity(problem_changed_metric2, @trustregions, Mchangedmetric2_handle, opts_adapt, Gmean_d, Amean_d);
        end
        info = add_time(info, time_precomp_A + time_precomp_G);
        results{5+idx_shift} = info;
        names{5+idx_shift} = name_prefix + rtr_type + " + Tools metric $\mathcal{P}^{(2)}X = K_0XG$" + name_postfix;
    end
end

%% Plot comparison
titletext = "Test Problem "+num2str(TP) + ":  m = " + num2str(m) + ",  n = " + num2str(n);
whatzoom_firstord = [0, ones(1,6)];
whatzoom_rtr = [0, ones(1, 4)];
what_zoom = [repmat(whatzoom_firstord, [1, length(cell_optsfirstord)]), ...
    repmat(whatzoom_rtr, [1, length(cell_optsrtr)])];
savename = "PAPER_stoch_diff_TP" + num2str(TP);
if rank_adaptivity
    fixed_rank = "rank_adaptive";
end
plot_manopt_results(results, names, fixed_rank, titletext, whatruns_final, what_zoom, savename)
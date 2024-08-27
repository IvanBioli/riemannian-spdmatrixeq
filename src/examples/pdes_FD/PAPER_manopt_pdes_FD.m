% Description: Finite difference discretization of 2D PDEs with separable variables
%   Executes experiments with Riemannian optimization

clearvars -except alpha n fixed_rank seed...
    sanity_checks whatruns_firstord whatruns_rtr ...
    rank_adaptivity opts_adaptivity...
    plot_and_save rel_tolgradnorm tolres debug; clc;
%% GENERATING PROBLEM INSTANCE
% Parameters of the problem
% alpha = 10;
n_terms = 3;
% n = 1000;

% Assembling the system
[data, u] = example1_paper(alpha, n_terms);
[A, B, C] = assemble_pdes_FD(n+1, data);

% General options
rhs_norm = stable_norm_fact(struct('U', C.L, 'V', C.R));
if exist("rel_tolgradnorm", "var")
    tolgradnorm = rel_tolgradnorm * rhs_norm;
else
    tolres = 1e-6 * rhs_norm;
    tolgradnorm = tolres / 2;
end
statsfun_ = @(problem, x, stats, store) sylv_statsfun(problem, x, stats, store);

% Options for conjugategradient
% RGD
opts_rgd = struct();
opts_rgd.beta_type = "steep";
opts_rgd.tolgradnorm = tolgradnorm;
opts_rgd.minstepsize = 1e-4 * opts_rgd.tolgradnorm;
opts_rgd.maxiter = 150;
opts_rgd.statsfun = statsfun_;
% RCG
opts_rcg = struct();
opts_rcg.beta_type = "H-S-SATO";
opts_rcg.tolgradnorm = tolgradnorm;
opts_rcg.minstepsize = 1e-4 * opts_rcg.tolgradnorm;
opts_rcg.maxiter = opts_rgd.maxiter;
opts_rcg.statsfun = statsfun_;
maxiter_kron = 100;
cell_optsfirstord = {opts_rgd, opts_rcg};

% Options for trustregions
opts_rtr_rhess = struct();
opts_rtr_rhess.projehess = false;
opts_rtr_rhess.tolgradnorm = tolgradnorm;
opts_rtr_rhess.maxiter = 150;
opts_rtr_rhess.maxinner = 250;
opts_rtr_rhess.statsfun = statsfun_;
opts_rtr_projehess = opts_rtr_rhess;
opts_rtr_projehess.projehess = true;
cell_optsrtr = {opts_rtr_rhess, opts_rtr_projehess};

% Set experiments to run
% sanity_checks = false;
% fixed_rank =  20;
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
    whatruns_firstord = ones(1, 14);
end
if ~exist("whatruns_rtr", "var")
    whatruns_rtr = ones(1, 5);
end
if n * fixed_rank >= 2000
    whatruns_firstord(2) = 0;
end
if n > 100
    whatruns_firstord([3, 9]) = 0;
end
whatruns_final = [];

%% RUN EXPERIMENTS
kmin = 1;
k = 8;

% Precomputations for ADI, P(X) = AXD + DXA preconditioner
[A_dom, B_dom, ~] = assemble_pdes_FD(n+1, data.dominant);
A_ = A_dom{3}; D_ = B_dom{3}; E_ = A_dom{4}; B_ = B_dom{4};
E_d = E_; D_d = D_;
time_shifts = tic();
b_dom = eigs(A_, E_, 1, 'largestabs', 'Tolerance', 1e-4);
a_dom = eigs(A_, E_, 1, 'smallestabs', 'Tolerance', 1e-4);
d_dom = eigs(B_, D_, 1, 'largestabs', 'Tolerance', 1e-4);
c_dom = eigs(B_, D_, 1, 'smallestabs', 'Tolerance', 1e-4);
time_shifts = toc(time_shifts);

% Precomputation for ADI, P(X) = TX + XT preconditioner
T = A_;
time_shifts_T = tic();
b = eigs(T, 1, 'largestabs', 'Tolerance', 1e-4);
a = eigs(T, 1, 'smallestabs', 'Tolerance', 1e-4);
time_shifts_T = toc(time_shifts_T);
% Below a second type of one-sided preconditioner, still not working well
% e = ones(n, 1); T = data.k_mean * spdiags([-e 2*e -e], -1:1, n, n);
% a = data.k_mean * 2 * (1 - cos(pi / (n+1))); % = eigs(T, 1, 'smallestreal');
% b = data.k_mean * 2 * (1 - cos(n * pi / (n+1))); % = eigs(T, 1, 'largestreal');

% Options for preconditioner of the form P(X) = AXD + DXA
AXDEXBprecond_opts = struct();
% Truncation options for Prec. Richardson and fADI
prec_rich_fADI_truncopts.method = ["relF", "absF"];
prec_rich_fADI_truncopts.tol = [0.1 * tolres / rhs_norm, 0.001 * tolres];


% DEFINE THE MANIFOLD, THE COST AND THE GRADIENT
% Problem with standard Euclidean metric
problem.A = A; problem.B = B; problem.C = C;
problem.M = fixedrankembeddedfactory(n, n, fixed_rank);
problem.cost = @(X, store) sylv_cost(X, A, B, C, store);
problem.grad = @(X, store) sylv_grad(X, A, B, C, store);
problem.rhs_norm = rhs_norm; problem.C = C;
problem.proj_ehess = @(X, eta, store) sylv_hess(X, eta, A, B, C, store, struct('projehess', true));
problem.linesearch = @(x, d, store) sylv_linesearch_init(x, d, store, problem);
problem_pcg = problem;
M_handle =@(r) fixedrankembeddedfactory(n, n, r);

% Problem with metric induced by P(X)=EXD
problem_changed_metric.A = A; problem_changed_metric.B = B; problem_changed_metric.C = C;
problem_changed_metric.M = fixedrank_metricEXD_factory(E_, D_, fixed_rank);
problem_changed_metric.cost = @(X, store) sylv_cost(X, A, B, C, store);
problem_changed_metric.grad = @(X, store) sylv_grad_EXD(X, A, B, C, store, E_d, D_d);
problem_changed_metric.rhs_norm = rhs_norm;
problem_changed_metric.proj_ehess = @(X, eta, store) sylv_hess_EXD(X, eta, A, B, C, store, E_d, D_d, struct('projehess', true));
problem_changed_metric.linesearch = @(x, d, store) sylv_linesearch_init(x, d, store, problem_changed_metric);
Mchangedmetric_handle =@(r) fixedrank_metricEXD_factory(E_, D_, r);

% Debug options
if exist("debug", "var")
    problem.debug = debug;
    problem_pcg.debug = debug;
    problem_changed_metric.debug = debug;
end

% Initialization point
init_norm = 1;
x0 = problem.M.rand(); x0.S = init_norm * x0.S;
x0_changed_metric = problem_changed_metric.M.SVD2repr(x0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RCG/RCG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:length(cell_optsfirstord)
    opts_firstord = cell_optsfirstord{j};
    opts_firstord_kron = opts_firstord;
    opts_firstord_kron.maxiter = maxiter_kron;
    idx_shift = (j - 1) * 14;
    % Get if RCG or RGD
    if opts_firstord.beta_type == "steep"|| opts_firstord.beta_type == "S-D"
        cg_type = "R-GD";
    else
        %cg_type = opts_firstord.beta_type + " R-NLCG";
        cg_type = " R-NLCG";
    end
    if length(whatruns_firstord) > 14
        whatruns_firstord_ = whatruns_firstord((j-1)*14+1 : j*14);
    else
        whatruns_firstord_ = whatruns_firstord;
    end
    whatruns_final = [whatruns_final, whatruns_firstord_];
    
    % R-GD/R-NLCG
    if whatruns_firstord_(1)
        if ~rank_adaptivity
            [~, ~, info, ~] = conjugategradient(problem, x0, opts_firstord);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = opts_firstord;
            [~, info] = sylv_rankadaptivity(problem, @conjugategradient, M_handle, opts_adapt);
        end
        results{1+idx_shift} = info;
        names{1+idx_shift} = name_prefix + cg_type + name_postfix;
    end
    
    % R-GD/R-NLCG + P(X) = A(X) (solution via vectorization)
    if whatruns_firstord_(2)
        problem_pcg.precon = @(X, H, store) sylv_precon_kron(X, H, A, B, C, store);
        
        if ~rank_adaptivity
            [~, ~, info, ~] = conjugategradient(problem_pcg, x0, opts_firstord_kron);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = opts_firstord_kron;
            [~, info] = sylv_rankadaptivity(problem_pcg, @conjugategradient, M_handle, opts_adapt);
        end
        results{2+idx_shift} = info;
        names{2+idx_shift} = name_prefix + cg_type + " + $\mathcal{P} = \mathcal{A}$ (kron)" + name_postfix;
    end
    
    % R-GD/R-NLCG + Preconditioned Richardson P(X) = TX + XT (exact Sylvester)
    if whatruns_firstord_(3)
        precond_options = struct();
        precond_options.sylv.A = full(T); precond_options.sylv.B = full(T);
        problem_pcg.precon = @(X, H, store) sylv_richardson_sylvesterexact(X, A, B, C, store, precond_options);
        if ~rank_adaptivity
            [~, ~, info, ~] = conjugategradient(problem_pcg, x0, opts_firstord);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = opts_firstord;
            [~, info] = sylv_rankadaptivity(problem_pcg, @conjugategradient, M_handle, opts_adapt);
        end
        results{3+idx_shift} = info;
        names{3+idx_shift} = name_prefix + cg_type + " + Prec. Rich. $\mathcal{P}^{(1)}X = AX+XB$ (Bartels Stewart)"+ name_postfix;
    end
    
    % Riemannian Truncated Preconditioned Richardson P(X) = TX + XT (kmin fADI step)
    if whatruns_firstord_(4)
        precond_options = struct();
        precond_options.fadi.A = T; precond_options.fadi.B = T;
        [p, q] = zolotarev_poles(kmin, a, b, a, b);
        precond_options.fadi.p = p; precond_options.fadi.q = q;
        precond_options.trunc_opts = prec_rich_fADI_truncopts;
        problem_pcg.precon = @(X, H, store) sylv_richardson_fADI(X, A, B, C, store, precond_options);
        
        if ~rank_adaptivity
            [~, ~, info, ~] = conjugategradient(problem_pcg, x0, opts_firstord);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = opts_firstord;
            [~, info] = sylv_rankadaptivity(problem_pcg, @conjugategradient, M_handle, opts_adapt);
        end
        info = add_time(info, time_shifts_T);
        results{4+idx_shift} = info;
        names{4+idx_shift} = name_prefix + cg_type + " + Prec. Rich. $\mathcal{P}^{(1)}X = AX+XB$ fADI$\times" +num2str(kmin) + "$"+ name_postfix;
    end
    
    % Riemannian Truncated Preconditioned Richardson P(X) = TX + XT (k fADI steps)
    if whatruns_firstord_(5)
        precond_options = struct();
        precond_options.fadi.A = T; precond_options.fadi.B = T;
        [p, q] = zolotarev_poles(k, a, b, a, b);
        precond_options.fadi.p = p; precond_options.fadi.q = q;
        precond_options.trunc_opts = prec_rich_fADI_truncopts;
        problem_pcg.precon = @(X, H, store) sylv_richardson_fADI(X, A, B, C, store, precond_options);
        
        if ~rank_adaptivity
            [~, ~, info, ~] = conjugategradient(problem_pcg, x0, opts_firstord);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = opts_firstord;
            [~, info] = sylv_rankadaptivity(problem_pcg, @conjugategradient, M_handle, opts_adapt);
        end
        info = add_time(info, time_shifts_T);
        results{5+idx_shift} = info;
        names{5+idx_shift} = name_prefix + cg_type + " + Prec. Rich. $\mathcal{P}^{(1)}X = AX+XB$ fADI$\times" +num2str(k) + "$"+ name_postfix;
    end
    
    % R-GD/R-NLCG + Riemannian preconditioner P(X) = TX + XT (exact)
    if whatruns_firstord_(6)
        problem_pcg.precon = @(X, H, store) sylv_precon_AXDEXB(X, H, T, T, store, AXDEXBprecond_opts);
        
        if ~rank_adaptivity
            [~, ~, info, ~] = conjugategradient(problem_pcg, x0, opts_firstord);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = opts_firstord;
            [~, info] = sylv_rankadaptivity(problem_pcg, @conjugategradient, M_handle, opts_adapt);
        end
        results{6+idx_shift} = info;
        names{6+idx_shift} = name_prefix + cg_type + " + $\mathcal{P}^{(1)}X = AX+XB$ (exact)" + name_postfix;
    end
    
    % R-GD/R-NLCG + Riemannian preconditioner P(X) = TX + XT (kmin tangADI steps)
    if whatruns_firstord_(7)
        precond_options = struct();
        [p, q] = zolotarev_poles(kmin, a, b, a, b);
        precond_options.p = p; precond_options.q = q;
        problem_pcg.precon = @(X, H, store) sylv_precon_tangADI(X, H, T, T, store, precond_options);
        
        if ~rank_adaptivity
            [~, ~, info, ~] = conjugategradient(problem_pcg, x0, opts_firstord);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = opts_firstord;
            [~, info] = sylv_rankadaptivity(problem_pcg, @conjugategradient, M_handle, opts_adapt);
        end
        info = add_time(info, time_shifts_T);
        results{7+idx_shift} = info;
        names{7+idx_shift} = name_prefix + cg_type + " + $\mathcal{P}^{(1)}X = AX+XB$ tangADI$\times" +num2str(kmin) + "$"+ name_postfix;
    end
    
    % R-GD/R-NLCG + Riemannian preconditioner P(X) = TX + XT (k tangADI steps)
    if whatruns_firstord_(8)
        precond_options = struct();
        [p, q] = zolotarev_poles(k, a, b, a, b);
        precond_options.p = p; precond_options.q = q;
        problem_pcg.precon = @(X, H, store) sylv_precon_tangADI(X, H, T, T, store, precond_options);
        
        if ~rank_adaptivity
            [~, ~, info, ~] = conjugategradient(problem_pcg, x0, opts_firstord);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = opts_firstord;
            [~, info] = sylv_rankadaptivity(problem_pcg, @conjugategradient, M_handle, opts_adapt);
        end
        info = add_time(info, time_shifts_T);
        results{8+idx_shift} = info;
        names{8+idx_shift} = name_prefix + cg_type + " + $\mathcal{P}^{(1)}X = AX+XB$ tangADI$\times" +num2str(k) + "$"+ name_postfix;
    end
    
    % Riemannian Truncated Preconditioned Richardson P(X) = AXD + DXA (exact Sylvester)
    if whatruns_firstord_(9)
        precond_options = struct();
        precond_options.sylv.A = full(A_); precond_options.sylv.D = full(D_);
        precond_options.sylv.E = full(E_);precond_options.sylv.B = full(B_);
        problem_pcg.precon = @(X, H, store) sylv_richardson_sylvesterexact(X, A, B, C, store, precond_options);
        
        if ~rank_adaptivity
            [~, ~, info, ~] = conjugategradient(problem_pcg, x0, opts_firstord);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = opts_firstord;
            [~, info] = sylv_rankadaptivity(problem_pcg, @conjugategradient, M_handle, opts_adapt);
        end
        results{9+idx_shift} = info;
        names{9+idx_shift} = name_prefix + cg_type + " + Prec. Rich. $\mathcal{P}^{(2)}X = AXD+EXB$ (Bartels Stewart)" + name_postfix;
    end
    
    % Riemannian Truncated Preconditioned Richardson P(X) = AXD + DXA (kmin fADI steps)
    if whatruns_firstord_(10)
        precond_options = struct();
        [p, q] = zolotarev_poles(kmin, a_dom, b_dom, c_dom, d_dom);
        precond_options.fadi.p = p; precond_options.fadi.q = q;
        precond_options.fadi.A = A_; precond_options.fadi.D = D_;
        precond_options.fadi.E = E_;precond_options.fadi.B = B_;
        precond_options.trunc_opts = prec_rich_fADI_truncopts;
        problem_pcg.precon = @(X, H, store) sylv_richardson_fADI(X, A, B, C, store, precond_options);
        
        if ~rank_adaptivity
            [~, ~, info, ~] = conjugategradient(problem_pcg, x0, opts_firstord);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = opts_firstord;
            [~, info] = sylv_rankadaptivity(problem_pcg, @conjugategradient, M_handle, opts_adapt);
        end
        info = add_time(info, time_shifts);
        results{10+idx_shift} = info;
        names{10+idx_shift} = name_prefix + cg_type + " + Prec. Rich. $\mathcal{P}^{(2)}X = AXD+EXB$ fADI$\times" +num2str(kmin) + "$"+ name_postfix;
    end
    
    % Riemannian Truncated Preconditioned Richardson P(X) = AXD + DXA (k fADI steps)
    if whatruns_firstord_(11)
        precond_options = struct();
        [p, q] = zolotarev_poles(k, a_dom, b_dom, c_dom, d_dom);
        precond_options.fadi.p = p; precond_options.fadi.q = q;
        precond_options.fadi.A = A_; precond_options.fadi.D = D_;
        precond_options.fadi.E = E_;precond_options.fadi.B = B_;
        precond_options.trunc_opts = prec_rich_fADI_truncopts;
        problem_pcg.precon = @(X, H, store) sylv_richardson_fADI(X, A, B, C, store, precond_options);
        
        if ~rank_adaptivity
            [~, ~, info, ~] = conjugategradient(problem_pcg, x0, opts_firstord);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = opts_firstord;
            [~, info] = sylv_rankadaptivity(problem_pcg, @conjugategradient, M_handle, opts_adapt);
        end
        info = add_time(info, time_shifts);
        results{11+idx_shift} = info;
        names{11+idx_shift} = name_prefix + cg_type + " + Prec. Rich. $\mathcal{P}^{(2)}X = AXD+EXB$ fADI$\times" +num2str(k) + "$"+ name_postfix;
    end
    
    % R-GD/R-NLCG + Riemannian preconditioner P(X) = AXD + DXA (exact solution via change of metric)
    if whatruns_firstord_(12)
        problem_changed_metric.precon = @(X, H, store) sylv_precon_AXDEXB(X, H, A_, B_, store, AXDEXBprecond_opts, D_, E_);
        opts_solver = opts_firstord;
        opts_solver.tolgradnorm = opts_solver.tolgradnorm / sqrt(normest(problem_changed_metric.M.E, 1e-4) * normest(problem_changed_metric.M.D, 1e-4));
        opts_solver.minstepsize = 1e-4 * opts_solver.tolgradnorm;
        if ~rank_adaptivity
            [~, ~, info, ~] = conjugategradient(problem_changed_metric, x0_changed_metric, opts_solver);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0_changed_metric;
            opts_adapt.opts_solver = opts_solver;
            [~, info] = sylv_rankadaptivity(problem_changed_metric, @conjugategradient, Mchangedmetric_handle, opts_adapt, D_d, E_d);
        end
        results{12+idx_shift} = info;
        names{12+idx_shift} = name_prefix + cg_type + " + $\mathcal{P}^{(2)}X = AXD+EXB$ (exact)" + name_postfix;
    end
    
    % R-GD/R-NLCG + Riemannian preconditioner P(X) = AXD + DXA (kmin tangADI steps)
    if whatruns_firstord_(13)
        precond_options = struct();
        [p, q] = zolotarev_poles(kmin, a_dom, b_dom, c_dom, d_dom);
        precond_options.p = p; precond_options.q = q;
        problem_pcg.precon = @(X, H, store) sylv_precon_tangADI(X, H, A_, B_, store, precond_options, E_, D_);
        
        if ~rank_adaptivity
            [~, ~, info, ~] = conjugategradient(problem_pcg, x0, opts_firstord);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = opts_firstord;
            [~, info] = sylv_rankadaptivity(problem_pcg, @conjugategradient, M_handle, opts_adapt);
        end
        info = add_time(info, time_shifts);
        results{13+idx_shift} = info;
        names{13+idx_shift} = name_prefix + cg_type + " + $\mathcal{P}^{(2)}X = AXD+EXB$ tangADI$\times" +num2str(kmin) + "$"+ name_postfix;
    end
    
    % R-GD/R-NLCG + Riemannian preconditioner P(X) = AXD + DXA (k tangADI steps)
    if whatruns_firstord_(14)
        precond_options = struct();
        [p, q] = zolotarev_poles(k, a_dom, b_dom, c_dom, d_dom);
        precond_options.p = p; precond_options.q = q;
        problem_pcg.precon = @(X, H, store) sylv_precon_tangADI(X, H, A_, B_, store, precond_options, E_, D_);
        
        if ~rank_adaptivity
            [~, ~, info, ~] = conjugategradient(problem_pcg, x0, opts_firstord);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = opts_firstord;
            [~, info] = sylv_rankadaptivity(problem_pcg, @conjugategradient, M_handle, opts_adapt);
        end
        info = add_time(info, time_shifts);
        results{14+idx_shift} = info;
        names{14+idx_shift} = name_prefix + cg_type + " + $\mathcal{P}^{(2)}X = AXD+EXB$ tangADI$\times" +num2str(k) + "$"+ name_postfix;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RTR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:length(cell_optsrtr)
    options_rtr = cell_optsrtr{j};
    idx_shift = length(cell_optsfirstord) * 14 + ...
        (j - 1) * 5;
    if length(whatruns_rtr) > 5
        whatruns_rtr_ = whatruns_rtr((j-1)*5+1 : j*5);
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
    problem_changed_metric.hess = @(X, eta, store) sylv_hess_EXD(X, eta, A, B, C, store, E_d, D_d, options_hess);
    
    % RTR
    if whatruns_rtr_(1)
        if ~rank_adaptivity
            [~, ~, info, ~] = trustregions(problem, x0, options_rtr);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = options_rtr;
            [~, info] = sylv_rankadaptivity(problem, @trustregions, M_handle, opts_adapt);
        end
        results{1+idx_shift} = info;
        names{1+idx_shift} = name_prefix + rtr_type + name_postfix;
    end
    
    % RTR + tCG + preconditioner P(X) = TX + XT (exact)
    if whatruns_rtr_(2)
        problem_pcg.precon = @(X, H, store) sylv_precon_AXDEXB(X, H, T, T, store, AXDEXBprecond_opts);
        
        if ~rank_adaptivity
            [~, ~, info, ~] = trustregions(problem_pcg, x0, options_rtr);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = options_rtr;
            [~, info] = sylv_rankadaptivity(problem_pcg, @trustregions, M_handle, opts_adapt);
        end
        results{2+idx_shift} = info;
        names{2+idx_shift} = name_prefix + rtr_type + " + $\mathcal{P}^{(1)}X = AX+XB$ (exact)" + name_postfix;
    end
    
    % RTR + tCG + preconditioner P(X) = TX + XT (k tangADI steps)
    if whatruns_rtr_(3)
        precond_options = struct();
        [p, q] = zolotarev_poles(k, a, b, a, b);
        precond_options.p = p; precond_options.q = q;
        problem_pcg.precon = @(X, H, store) sylv_precon_tangADI(X, H, T, T, store, precond_options);
        
        if ~rank_adaptivity
            [~, ~, info, ~] = trustregions(problem_pcg, x0, options_rtr);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = options_rtr;
            [~, info] = sylv_rankadaptivity(problem_pcg, @trustregions, M_handle, opts_adapt);
        end
        info = add_time(info, time_shifts_T);
        results{3+idx_shift} = info;
        names{3+idx_shift} = name_prefix + rtr_type + " + $\mathcal{P}^{(1)}X = AX+XB$ tangADI$\times" +num2str(k) + "$"+ name_postfix;
    end
    
    % RTR + tCG + preconditioner P(X) = AXD + DXA (exact solution via change of metric)
    if whatruns_rtr_(4)
        problem_changed_metric.precon = @(X, H, store) sylv_precon_AXDEXB(X, H, A_, B_, store, AXDEXBprecond_opts, D_, E_);
        opts_solver = options_rtr;
        opts_solver.tolgradnorm = opts_solver.tolgradnorm / sqrt(normest(problem_changed_metric.M.E, 1e-4) * normest(problem_changed_metric.M.D, 1e-4));
        opts_solver.minstepsize = 1e-4 * opts_solver.tolgradnorm;
        if ~rank_adaptivity
            [~, ~, info, ~] = trustregions(problem_changed_metric, x0_changed_metric, opts_solver);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0_changed_metric;
            opts_adapt.opts_solver = opts_solver;
            [~, info] = sylv_rankadaptivity(problem_changed_metric, @trustregions, Mchangedmetric_handle, opts_adapt, D_d, E_d);
        end
        results{4+idx_shift} = info;
        names{4+idx_shift} = name_prefix + rtr_type + " + $\mathcal{P}^{(2)}X = AXD+EXB$ (exact)" + name_postfix;
    end
    
    % RTR + tCG + preconditioner P(X) = AXD + DXA (k tangADI steps)
    if whatruns_rtr_(5)
        precond_options = struct();
        [p, q] = zolotarev_poles(k, a_dom, b_dom, c_dom, d_dom);
        precond_options.p = p; precond_options.q = q;
        problem_pcg.precon = @(X, H, store) sylv_precon_tangADI(X, H, A_, B_, store, precond_options, E_, D_);
        
        if ~rank_adaptivity
            [~, ~, info, ~] = trustregions(problem_pcg, x0, options_rtr);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = options_rtr;
            [~, info] = sylv_rankadaptivity(problem_pcg, @trustregions, M_handle, opts_adapt);
        end
        info = add_time(info, time_shifts_T);
        results{5+idx_shift} = info;
        names{5+idx_shift} = name_prefix + rtr_type + " + $\mathcal{P}^{(2)}X = AXD+EXB$ tangADI$\times" +num2str(k) + "$"+ name_postfix;
    end
end

%%
if sanity_checks
    whatruns_final = [whatruns_final, 1, 1, 1];
    
    % DEFINE THE MANIFOLD, THE COST AND THE GRADIENT
    problem_check = struct();
    problem_check.M = euclideanfactory(n,n);
    problem_check.costgrad = @(X) gen_sylv_posdef_fr_costgrad(X, A, B, C.L * C.R');
    problem_check.rhs_norm = problem.rhs_norm;
    X0 = x0.U * x0.S * x0.V';
    
    % Steepest descent
    [~, ~, info, ~] = conjugategradient(problem_check, X0, opts_rgd);
    results{length(whatruns)+1} = info;
    names{length(whatruns)+1} = "Full rank GD";
    
    % Preconditioned steepest descent
    problem_check.precon =@(X,H) sylvester(full(T), full(T), H);
    [~, ~, info, ~] = conjugategradient(problem_check, X0, opts_rgd);
    results{length(whatruns)+2} = info;
    names{length(whatruns)+2} = "Full rank pGD $\mathcal{P}^{(1)}X = AX+XB$";
    
    % Preconditioned steepest descent
    problem_check.precon =@(X,H) bartelsStewart(A_, D_, E_, B_, H, false, false);
    [~, ~, info, ~] = conjugategradient(problem_check, X0, opts_rgd);
    results{length(whatruns)+3} = info;
    names{length(whatruns)+3} = "Full rank pGD $\mathcal{P}^{(2)}X = AXD+EXB$";
end
%% Plot comparison
if ~exist("plot_and_save", "var")
    plot_and_save = true;
end
if plot_and_save
    titletext = "$a =" + alpha + ", n=" + n + ", r = " + fixed_rank + "$";
    whatzoom_firstord = ones(1, 14); whatzoom_firstord([4, 7, 10, 13]) = 0;
    whatzoom_rtr = [0, ones(1,4)];
    whatzoom = [repmat(whatzoom_firstord, [1, length(cell_optsfirstord)]), ...
        repmat(whatzoom_rtr, [1, length(cell_optsrtr)])];
    if sanity_checks
        whatzoom = [whatzoom, 0, 1, 1];
    end
    savename = "PAPER_pdes_FD_n" + num2str(n) + "_a" + num2str(alpha);
    if rank_adaptivity
        fixed_rank = "rank_adaptive";
    end
    plot_manopt_results(results, names, fixed_rank, titletext, whatruns_final, whatzoom, savename)
end
% Description: Reachability Gramian of a bilinear control system
%   Executes experiments with Riemannian optimization

clearvars -except seed k_rail fixed_rank ...
    whatruns_firstord whatruns_rtr ...
    rank_adaptivity opts_adaptivity...
    plot_and_save rel_tolgradnorm tolres debug k_tangADI; clc;
%% GENERATING PROBLEM INSTANCE
% Parameters of the problem
% k_rail = 5; fixed_rank = 245;
load_rail;

% General options
rhs_norm = stable_norm_fact(struct('V', B, 'D', eye(size(B,2))));
if exist("rel_tolgradnorm", "var")
    tolgradnorm = rel_tolgradnorm * rhs_norm;
else
    tolres = 1e-6 * rhs_norm;
    tolgradnorm = tolres / 1.5;
end
statsfun_ = @(problem, x, stats, store) lyap_psd_statsfun(problem, x, stats, store);

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
    whatruns_firstord = [zeros(1,7), zeros(1, 7)];
    whatruns_firstord([14]) = 1;
end
if ~exist("whatruns_rtr", "var")
    whatruns_rtr = [zeros(1, 5), zeros(1,5)];
    % whatruns_rtr([9]) = 1;
end
whatruns_final = [];
%% RUN EXPERIMENTS
k_fadi = 8;
if ~exist("k_tangADI", "var")
    k_tangADI = 4;
end

% Precomputations for the ADI shifts
I = speye(n);
time_shifts_onesided = tic();
a_onesided = eigs(A, 1, 'smallestabs');
b_onesided = eigs(A, 1, 'largestabs');
time_shifts_onesided = toc(time_shifts_onesided);

time_shifts_twosided = tic();
a_twosided = eigs(A, M, 1, 'smallestabs');
b_twosided = eigs(A, M, 1, 'largestabs');
time_shifts_twosided = toc(time_shifts_twosided);

% Options for preconditioner of the form P(X) = AXM + MXA
AXMMXAprecond_opts = struct();
AXMMXAprecond_opts.small_memory = true;
% AXMMXAprecond_opts.assemble_F = false;
% Truncation options for Prec. Richardson and fADI
prec_rich_fADI_truncopts.method = ["relF", "absF"];
prec_rich_fADI_truncopts.tol = [0.1 * tolres / rhs_norm, 0.0001 * tolres];

% DEFINE THE MANIFOLD, THE COST AND THE GRADIENT
% Problem with standard Euclidean metric
problem.A = A; problem.Mass = M; problem.B = B; problem.N = N;
problem.M = sympsdfixedrankembeddedfactory(n, fixed_rank);
problem.cost = @(X, store) lyap_psd_cost(X, A, M, N, B, store);
problem.grad = @(X, store) lyap_psd_grad(X, A, M, N, B, store);
problem.rhs_norm = rhs_norm;
problem.proj_ehess = @(X, eta, store) lyap_psd_hess(X, eta, A, M, N, B, store, struct('projehess', true));
problem.linesearch = @(x, d, store) lyap_psd_linesearch_init(x, d, store, problem);
problem_pcg = problem;
M_handle =@(r) sympsdfixedrankembeddedfactory(n, r);

% Problem with metric induced by P(X)=MXM
problem_changed_metric.A = A; problem_changed_metric.Mass = M; problem_changed_metric.B = B; problem_changed_metric.N = N;
problem_changed_metric.M = sympsdfixedrank_metricMXM_factory(M, fixed_rank);
problem_changed_metric.cost = @(X, store) lyap_psd_cost(X, A, M, N, B, store);
problem_changed_metric.grad = @(X, store) lyap_psd_grad_MXM(X, A, M, N, B, store, M_d);
problem_changed_metric.rhs_norm = rhs_norm;
problem_changed_metric.proj_ehess = @(X, eta, store) lyap_psd_hess_MXM(X, eta, A, M, N, B, store, struct('projehess', true), M_d);
problem_changed_metric.linesearch = @(x, d, store) lyap_psd_linesearch_init(x, d, store, problem_changed_metric);
Mchangedmetric_handle =@(r) sympsdfixedrank_metricMXM_factory(M, r);

% Do some steps of fADI as initialization
start_init = tic();
p = zolotarev_poles(ceil(fixed_rank / size(B,2)), a_twosided, b_twosided);
Z = fadi_lyap(A, B, p, M);
[Q,R] = qr(Z, 0); [V,D] = eig(symmetrize(R*R')); V = Q*V;
[~,ind] = sort(diag(D), 'descend'); V = V(:, ind); D = D(ind, ind);
x0.V = V(:,1:fixed_rank); x0.D = D(1:fixed_rank,1:fixed_rank);
time_init = toc(start_init);
x0_changed_metric = problem_changed_metric.M.eig2repr(x0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RCG/RCG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:length(cell_optsfirstord)
    opts_firstord = cell_optsfirstord{j};
    opts_firstord_kron = opts_firstord;
    opts_firstord_kron.maxiter = maxiter_kron;
    idx_shift = (j - 1) * 7;
    % Get if RCG or RGD
    if opts_firstord.beta_type == "steep"|| opts_firstord.beta_type == "S-D"
        cg_type = "R-GD";
    else
        %cg_type = opts_firstord.beta_type + " R-NLCG";
        cg_type = " R-NLCG";
    end
    if length(whatruns_firstord) > 7
        whatruns_firstord_ = whatruns_firstord((j-1)*7+1 : j*7);
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
            [~, info] = lyap_psd_rankadaptivity(problem, @conjugategradient, M_handle, opts_adapt);
        end
        info = add_time(info, time_init + time_shifts_twosided);
        results{1+idx_shift} = info;
        names{1+idx_shift} = name_prefix + cg_type + name_postfix;
    end
    
    % Riemannian Truncated Preconditioned Richardson P(X) = AX + XA (k_fadi fADI steps)
    if whatruns_firstord_(2)
        precond_options = struct();
        precond_options.fadi.A = A;
        precond_options.fadi.p = zolotarev_poles(k_fadi, a_onesided, b_onesided);
        precond_options.trunc_opts = prec_rich_fADI_truncopts;
        problem_pcg.precon = @(X, eta, store) lyap_psd_richardson_ADI(X, A, M, N, B, store, precond_options);
        
        if ~rank_adaptivity
            [~, ~, info, ~] = conjugategradient(problem_pcg, x0, opts_firstord);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = opts_firstord;
            [~, info] = lyap_psd_rankadaptivity(problem_pcg, @conjugategradient, M_handle, opts_adapt);
        end
        results{2+idx_shift} = info;
        info = add_time(info, time_init);
        info = add_time(info, time_shifts_onesided);
        names{2+idx_shift} = name_prefix + cg_type + " + Prec. Rich. $\mathcal{P}^{(1)}X = AX+XA$ fADI$\times" +num2str(k_fadi) + "$"+ name_postfix;
    end
    
    % R-GD/R-NLCG + Riemannian preconditioner P(X) = AX + XA (exact)
    if whatruns_firstord_(3)
        problem_pcg.precon = @(X, eta, store) lyap_psd_precon_AXMMXA(X, eta, A, store, AXMMXAprecond_opts);
        
        if ~rank_adaptivity
            [~, ~, info, ~] = conjugategradient(problem_pcg, x0, opts_firstord);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = opts_firstord;
            [~, info] = lyap_psd_rankadaptivity(problem_pcg, @conjugategradient, M_handle, opts_adapt);
        end
        info = add_time(info, time_init + time_shifts_twosided);
        results{3+idx_shift} = info;
        names{3+idx_shift} = name_prefix + cg_type + " + $\mathcal{P}^{(1)}X = AX+XA$ (exact)" + name_postfix;
    end
    
    % R-GD/R-NLCG + Riemannian preconditioner P(X) = AX + XA (k_tangADI tangADI steps)
    if whatruns_firstord_(4)
        precond_options = struct();
        precond_options.p = zolotarev_poles(k_tangADI, a_onesided, b_onesided);
        problem_pcg.precon = @(X, eta, store) lyap_psd_precon_tangADI(X, eta, A, store, precond_options);
        
        if ~rank_adaptivity
            [~, ~, info, ~] = conjugategradient(problem_pcg, x0, opts_firstord);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = opts_firstord;
            [~, info] = lyap_psd_rankadaptivity(problem_pcg, @conjugategradient, M_handle, opts_adapt);
        end
        results{4+idx_shift} = info;
        info = add_time(info, time_init);
        info = add_time(info, time_shifts_onesided);
        names{4+idx_shift} = name_prefix + cg_type + " + $\mathcal{P}^{(1)}X = AX+XA$ tangADI$\times" +num2str(k_tangADI) + "$"+ name_postfix;
    end
    
    % Riemannian Truncated Preconditioned Richardson P(X) = AXM + MXA (k_faDI fADI steps)
    if whatruns_firstord_(5)
        precond_options = struct();
        precond_options.fadi.A = A;
        precond_options.fadi.M = M;
        precond_options.fadi.p = zolotarev_poles(k_fadi, a_twosided, b_twosided);
        precond_options.trunc_opts = prec_rich_fADI_truncopts;
        problem_pcg.precon = @(X, eta, store) lyap_psd_richardson_ADI(X, A, M, N, B, store, precond_options);
        
        if ~rank_adaptivity
            [~, ~, info, ~] = conjugategradient(problem_pcg, x0, opts_firstord);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = opts_firstord;
            [~, info] = lyap_psd_rankadaptivity(problem_pcg, @conjugategradient, M_handle, opts_adapt);
        end
        info = add_time(info, time_init);
        info = add_time(info, time_shifts_twosided);
        results{5+idx_shift} = info;
        names{5+idx_shift} = name_prefix + cg_type + " + Prec. Rich. $\mathcal{P}^{(2)}X = AXM+MXA$ fADI$\times" +num2str(k_fadi) + "$"+ name_postfix;
    end
    
    % R-GD/R-NLCG + Riemannian preconditioner P(X) = AXM + MXA (exact solution via change of metric)
    if whatruns_firstord_(6)
        problem_changed_metric.precon = @(X, eta, store) lyap_psd_precon_AXMMXA(X, eta, A, store, AXMMXAprecond_opts, M);
        opts_solver = opts_firstord;
        opts_solver.tolgradnorm = opts_solver.tolgradnorm / normest(problem_changed_metric.M.Mass, 1e-4);
        opts_solver.minstepsize = 1e-4 * opts_solver.tolgradnorm;
        if ~rank_adaptivity
            [~, ~, info, ~] = conjugategradient(problem_changed_metric, x0_changed_metric, opts_solver);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0_changed_metric;
            opts_adapt.opts_solver = opts_solver;
            [~, info] = lyap_psd_rankadaptivity(problem_changed_metric, @conjugategradient, Mchangedmetric_handle, opts_adapt, D_d, E_d);
        end
        info = add_time(info, time_init + time_shifts_twosided);
        results{6+idx_shift} = info;
        names{6+idx_shift} = name_prefix + cg_type + " + $\mathcal{P}^{(2)}X = AXM+MXA$ (exact)" + name_postfix;
    end
    
    % R-GD/R-NLCG + Riemannian preconditioner P(X) = AXM + MXA (k_tangADI tangADI steps)
    if whatruns_firstord_(7)
        precond_options = struct();
        precond_options.p = zolotarev_poles(k_tangADI, a_twosided, b_twosided);
        problem_pcg.precon = @(X, eta, store) lyap_psd_precon_tangADI(X, eta, A, store, precond_options, M);
        
        if ~rank_adaptivity
            [~, ~, info, ~] = conjugategradient(problem_pcg, x0, opts_firstord);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = opts_firstord;
            [~, info] = lyap_psd_rankadaptivity(problem_pcg, @conjugategradient, M_handle, opts_adapt);
        end
        info = add_time(info, time_init);
        info = add_time(info, time_shifts_twosided);
        results{7+idx_shift} = info;
        names{7+idx_shift} = name_prefix + cg_type + " + $\mathcal{P}^{(2)}X = AXM+MXA$ tangADI$\times" +num2str(k_tangADI) + "$"+ name_postfix;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RTR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:length(cell_optsrtr)
    options_rtr = cell_optsrtr{j};
    idx_shift = length(cell_optsfirstord) * 7 + ...
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
    problem.hess = @(X, eta, store) lyap_psd_hess(X, eta, A, M, N, B, store, options_hess);
    problem_pcg = problem;
    problem_changed_metric.hess = @(X, eta, store) lyap_psd_hess_MXM(X, eta, A, M, N, B, store, options_hess, M_d);
    
    % RTR
    if whatruns_rtr_(1)
        if ~rank_adaptivity
            [~, ~, info, ~] = trustregions(problem, x0, options_rtr);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = options_rtr;
            [~, info] = lyap_psd_rankadaptivity(problem, @trustregions, M_handle, opts_adapt);
        end
        info = add_time(info, time_init + time_shifts_twosided);
        results{1+idx_shift} = info;
        names{1+idx_shift} = name_prefix + rtr_type + name_postfix;
    end
    
    % RTR + tCG + preconditioner P(X) = AX + XA (exact)
    if whatruns_rtr_(2)
        problem_pcg.precon = @(X, eta, store) lyap_psd_precon_AXMMXA(X, eta, A, store, AXMMXAprecond_opts);
        
        if ~rank_adaptivity
            [~, ~, info, ~] = trustregions(problem_pcg, x0, options_rtr);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = options_rtr;
            [~, info] = lyap_psd_rankadaptivity(problem_pcg, @trustregions, M_handle, opts_adapt);
        end
        info = add_time(info, time_init + time_shifts_twosided);
        results{2+idx_shift} = info;
        names{2+idx_shift} = name_prefix + rtr_type + " + $\mathcal{P}^{(1)}X = AX+XA$ (exact)" + name_postfix;
    end
    
    % RTR + tCG + preconditioner P(X) = AX + XA (k_tangADI tangADI steps)
    if whatruns_rtr_(3)
        precond_options = struct();
        precond_options.p = zolotarev_poles(k_tangADI, a_onesided, b_onesided);
        problem_pcg.precon = @(X, eta, store) lyap_psd_precon_tangADI(X, eta, A, store, precond_options);
        
        if ~rank_adaptivity
            [~, ~, info, ~] = trustregions(problem_pcg, x0, options_rtr);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = options_rtr;
            [~, info] = lyap_psd_rankadaptivity(problem_pcg, @trustregions, M_handle, opts_adapt);
        end
        info = add_time(info, time_init + time_shifts_twosided);
        results{3+idx_shift} = info;
        names{3+idx_shift} = name_prefix + rtr_type + " + $\mathcal{P}^{(1)}X = AX+XA$ (exact)" + name_postfix;
    end
    
    % RTR + tCG + preconditioner P(X) = AXM + MXA (exact solution via change of metric)
    if whatruns_rtr_(4)
        problem_changed_metric.precon = @(X, eta, store) lyap_psd_precon_AXMMXA(X, eta, A, store, AXMMXAprecond_opts, M);
        opts_solver = opts_firstord;
        opts_solver.tolgradnorm = opts_solver.tolgradnorm / normest(problem_changed_metric.M.Mass, 1e-4);
        opts_solver.minstepsize = 1e-4 * opts_solver.tolgradnorm;
        
        if ~rank_adaptivity
            [~, ~, info, ~] = trustregions(problem_changed_metric, x0_changed_metric, opts_solver);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0_changed_metric;
            opts_adapt.opts_solver = opts_solver;
            [~, info] = lyap_psd_rankadaptivity(problem_changed_metric, @trustregions, Mchangedmetric_handle, opts_adapt, D_d, E_d);
        end
        info = add_time(info, time_init + time_shifts_twosided);
        results{4+idx_shift} = info;
        names{4+idx_shift} = name_prefix + rtr_type + " + $\mathcal{P}^{(2)}X = AXM+MXA$ (exact)" + name_postfix;
    end
    
    % RTR + tCG + preconditioner P(X) = AXM + MXA (k_tangADI tangADI steps)
    if whatruns_rtr_(5)
        precond_options = struct();
        precond_options.p = zolotarev_poles(k_tangADI, a_twosided, b_twosided);
        problem_pcg.precon = @(X, eta, store) lyap_psd_precon_tangADI(X, eta, A, store, precond_options, M);
        
        if ~rank_adaptivity
            [~, ~, info, ~] = trustregions(problem_pcg, x0, options_rtr);
        else
            opts_adapt = opts_adaptivity;
            opts_adapt.X0 = x0;
            opts_adapt.opts_solver = options_rtr;
            [~, info] = lyap_psd_rankadaptivity(problem_pcg, @trustregions, M_handle, opts_adapt);
        end
        info = add_time(info, time_init + time_shifts_twosided);
        results{5+idx_shift} = info;
        names{5+idx_shift} = name_prefix + rtr_type + " + $\mathcal{P}^{(2)}X = AXM+MXA$ tangADI$\times" +num2str(k_tangADI) + "$"+ name_postfix;
    end
end

%% Plot comparison
if ~exist("plot_and_save", "var")
    plot_and_save = true;
end
if plot_and_save
    titletext = "Rail profile problem, $n = " + n + "$";
    savename = "rail_n" + num2str(n) + "_tangADIx" + num2str(k_tangADI);
    whatzoom = [];
    if rank_adaptivity
        fixed_rank = "rank_adaptive";
    end
    plot_manopt_results(results, names, fixed_rank, titletext, whatruns_final, whatzoom, savename)
end
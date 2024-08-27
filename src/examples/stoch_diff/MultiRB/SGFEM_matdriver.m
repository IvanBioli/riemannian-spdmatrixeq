%  ---------- Driver for SGFEM matrix equations --------------------------
%
% Performs essential pre-processing of SGFEM matrices + vectors. Selects 
% parameters (s_parameter) + other options for MultiRB solver and then calls
% MultiRB to solve matrix equation. 
%
% 'Pre-processing' involves preconditioning (here, based on mean-stiffness
% matrix K_0) + selecting a set of shifts (alphas). 
%
% NB: Shifts and parameters are problem-dependent and not optimized here. 
% ***This driver should ONLY be used for the test problems provided***. 
%
% -------------------------------------------------------------------------
% REFERENCE: An Efficient Reduced Basis Solver for Stochastic Galerkin
% Matrix Equations, C.E. Powell, D. Silvester, V. Simoncini, 
% SIAM J. Sci. Comput., Vol 39, No. 1, pp. A141--A163, (2017). 
% -------------------------------------------------------------------------
%
% Data for two test problems is provided: TP_five.mat (Example 5.1) and 
% TP_two.mat (Example 5.2). To run driver, just type
%
% >> SGFEM_matdriver
% 
% User is prompted to choose a test problem and then a method for selecting
% the parameters (s_parameter) for MultiRB.
% 
% >> Parameter-free or parameter-dependent version (1/2):
%
% (See paper for details). No default set - you must specify an option! 
% Parameters are chosen differently for each test problem.
%
% -------------------------------------------------------------------------
% SGFEM discretisation: Q1 FEM on uniform grid (grid level ell) + global 
% polynomials of total degree p or less in M variables. PDE is elliptic so
% these matrix equations correspond to linear systems Ax=b with A SPD. All
% G and K matrices are SYMMETRIC and K_0 is positive definite.   
% -------------------------------------------------------------------------
% Test Problem 5 (TP_five.mat, ell = 7, p=5, m =9) 
% Dim K matrices: 16,129, Dim G matrices = 2,002, 
% No. terms in matrix equation (m+1): 10  
% Total no. equations (Ax=b): 32,290,258.
% --------------------------
% Test Problem 2 (TP_two.mat, ell = 7, p=5, m = 8, sigma=0.3)  
% Dim_K matrices: 16,129, Dim G matrices: 1,287, 
% No. terms in matrix equation (m+1): 9  
% Total no. equations (Ax=b): 20,758,023.
% ---------------------------
% Code to generate other SGFEM elliptic test problems available at:
% S-IFISS toolbox:    www.manchester.ac.uk/ifiss/s-ifiss.html
% ------------------------------------------------------------------------
%
% Copyright (c): C.E. Powell, V. Simoncini, 12th August 2019.
%
% ------------------------------------------------------------------------

%TP=input('\nTest problem 5 (Example 5.1) or 2 (Example 5.2)? (5/2):');
if TP==5
    load TP_five
elseif TP==2
    load TP_two
end

fprintf('\n Dim_K = %d,  Dim_G = %d, no. terms = %d  \n', size(K{1},1),size(G{1},1),length(K))
fprintf('\n No. equations (Ax=b): %d  \n', size(K{1},1)*size(G{1},1))

n=size(K{1},1); m=size(G{1},1);
start = tic();
% Permute/sort mean stiffness matrices, based on Amean.
psort=symamd(Amean);
Amean_p=Amean(psort,psort);
for ind=1:length(K)
    K_p{ind}=K{ind}(psort,psort);
end

rhs1=fnew(1:n); rhs2=eye(m,1); rhs1_p=rhs1(psort); 
% Cholesky factorisation of sorted Amean
R=chol(Amean_p,'lower'); 
I=speye(n); Gnew=G;

% ----------------Choose shifts ------------------------------------
opts.tol=1e-4;  % NB: 'eigs' is VERY sensitive to this.
% Make sure K-matrices are symmetric (they should be).
K_p{1}=(K_p{1}+K_p{1}')/2;
K_p{2}=(K_p{2}+K_p{2}')/2;

% This could be done more cheaply with data from a COARSER problem
% min and max eigs of hat{K}_{1} - 1st term could be strictly positive
emin2=eigs(K_p{2},K_p{1},1,'SA',opts);
emax2=eigs(K_p{2},K_p{1},1,'LA',opts);
n_m=size(K_p,2);
emean=mean([(emin2),abs(emax2)]);
alphas(1)=1-emean;
alphas(2:n_m)=(1)*ones(n_m-1,1);        % center all other spectra around 1
    
%shifts = alphas'                       % Uncomment to print to screen

% Modify the orginal (sorted) Ks and Gs (apply shifts) 
 for k=2:size(K_p,2)
     K_p{k}=(K_p{k}+K_p{k}')/2; % Make sure symmetric
     K_p{k}=K_p{k}+alphas(k-1)*K_p{1};
     Gnew{1}=Gnew{1}-alphas(k-1)*Gnew{k};
 end
 
% --------------Choose parameters s for MultiRB -------------------

% p_option=input('\n Parameter-free or parameter-dependent version (1/2):');
if p_option==1
    % This makes it equivalent to 'parameter-free version' for above alphas;
    s_parameter(1)=2-alphas(1);
    s_parameter(2)=1;
    %s_parameter                             % Uncomment to print to screen
elseif p_option==2
    %TP=input('\nTest problem 2 or 5? (2/5):');
    if TP==5
        % For Test Problem 5 only (K_1 is indefinite)
        aa=emin2+alphas(1); bb=emax2+alphas(1);   
        S_interval=[aa,bb];
        s_nodes=5;                           % Choose 5 nodes (could vary)
        snew=get_nodes2(aa,bb,s_nodes);      % Use interval for A_1;
        s_parameter=snew;                    % Uncomment to print to screen
    elseif TP==2
        % For Test Problem 2 only (K_1 is positive definite)
        % Get min and max eigs of hat{K}_{2} 
        emin3=eigs(K_p{3},K_p{1},1,'SA',opts);
        emax3=eigs(K_p{3},K_p{1},1,'LA',opts);
        S_interval=[emin3,emax3];
        s_nodes=5;                             % Choose 5 nodes (could vary)
        snew=get_nodes2(emin3,emax3,s_nodes);  % Use interval for A_2;
        s_parameter=snew;                      % Uncomment to print to screen
    end
end
time_precomp = toc(start);

% ------------ Options for MultiRB --------------------------------
res_method='4';                  % can't change this
rat_solve='1'; mmax=200;         % can change
param.max_space_dim=mmax;
param.period=1;                  % can change
param.rat_solve=rat_solve;
param.res_method=res_method;
tic;
fprintf('\n -------------------- MultiRB Solve ------------------------\n \n')

[X1,X2,dimV,final_err,avg_inner,error_vec,iv_vec, info]=MultiRB(K_p,Gnew,rhs1_p,rhs2,Amean_p,R,param,s_parameter,tol);
etoc=toc; 
info = add_time(info, time_precomp);
X1(psort,:)=X1;
fprintf('\n Total execution time: %9.4e seconds \n',etoc)
fprintf('\n ----------------------------------------------------------\n \n')
fprintf('no_terms   dim_K   dim_G   n_k   Rank    final_err   avg_inner  time(s)  \n')
fprintf('\n  %2d       %3d   %d    %3d    %2d    %9.4e    %4.2f     %9.4e  \n \n', [noarv+1, n,m, dimV, size(X1,2), final_err   avg_inner  etoc])
fprintf('----------------------------------------------------------\n \n')
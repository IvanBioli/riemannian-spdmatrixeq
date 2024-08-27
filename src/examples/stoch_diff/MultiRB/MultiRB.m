function [X1,X2,dimV,final_err,avg_inner,error_vec,iv_vec,info]=MultiRB(K,G,rhs1,rhs2,M,M1,param,snew,tol)
%
% ---------------------- MultiRB solver -----------------------------------
%
% function [X1,X2,dimV,final_err,avg_inner,error_vec,iv_vec]=MultiRB(K,G,rhs1,rhs2,M,M1,param,snew);
%
% Solve  A x = b  arising from SGFEM stochastic diffusion problem using 
% matrix equation formulation:
%
%         K0 X G0 + K1 X G1 + .... + Km X Gm = rhs1*rhs2'
%
% Note - problem is pre-processed (preconditioned + shifted) first by
% running SGFEM_matriver.m
%
% INPUTS:
% K            set of m+1 stiffness matrices
% G            set of m+1 stochastic matrices
% rhs1,rhs2    such that  b=vec(rhs1*rhs2')
% M, M1        preconditioner, M = M1 * M1'
% param        contains MultiRB options, e,g:
%              param.max_space_dim=200;  (max approxn space dimension)
%              param.period=1;  (how often to check residual (no. its))
%              param.res_method='4'; (stopping criterion:  4 = relative diff soln)
%              param.rat_solve='1';  (RKS solver: 0=direct,1=iterative)
%
% OUTPUTS:
% X1, X2        X = X1*X2' final approximate solution
% dimV          dimension of final approximation space
% final_err     final relative difference in solution (see param.res_method)
% avg_inner     average no. inner iterations to solve projected problem
%
% -------------------------------------------------------------------------
% REFERENCE: An Efficient Reduced Basis Solver for Stochastic Galerkin
% Matrix Equations, C.E. Powell, D. Silvester, V. Simoncini, 
% SIAM J. Sci. Comput., Vol 39, No. 1, pp. A141--A163, (2017). 
% -------------------------------------------------------------------------
% 
% Copyright (c): V. Simoncini, C.E. Powell, 12th August 2019.
%
% ------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%% ADDED %%%%%%%%%%%%%%%%%%%%%%%%%%%
timetic = tic;
C.L = rhs1; C.R = rhs2;
rhs_norm = stable_norm_fact(struct('U', rhs1, 'V', rhs2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

compute_period=param.period;
rat_solve=param.rat_solve;
res_method=param.res_method;
mmax=param.max_space_dim;
nterm=size(G,2);
 

tol_drop=.99;  % controls how many basis vectors to add at each iteration
nofirst=0;
Y=0;
Mrhs1=M1\rhs1;
normMrhs1=norm(Mrhs1);
V = (Mrhs1)/normMrhs1;   % Include rhs in the basis.

tot_inner=[];
W=speye(size(G{1}));
rhs1m=V'*Mrhs1;
rhs2m=W'*rhs2;
normR0=norm(Mrhs1)*norm(rhs2);
normR0_noprec=norm(rhs1)*norm(rhs2);
nrmres=1;
nrmres_noprec=1;
error_vec(1)=nrmres_noprec;

Y0=[];
Gm=G;

s_nodes=length(snew);  % number of 's' parameters. 

rat_space=1;
I=speye(size(K{1},1));
wrk=M1'\V;
% Projected matrix
for ind=1:nterm
    Km{ind}=wrk'*(K{ind}*wrk);  %(M1'\V)));
end
m=size(G{1},1); n=size(K{1},1);
i=0; tot_it=0; ir=0; Y=[];

fprintf('      no.its      n_k       G_dim    error (rel diff)   no.inner its\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%% ADDED %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save stats in a struct array info, and preallocate.
stats = savestats();
info(1) = stats;
info(min(10000, mmax+1)).iter = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% main iteration
while (i < mmax & nrmres_noprec>tol) 
    i=i+1;
    tot_it=tot_it+1;
    if (i > size(V,2)), fprintf('exausted approx space\n'),end
    % Get new basis vectors
    if (rat_space)                % rational space
         ir=ir+1; if (ir>s_nodes),ir=1;end
         if (rat_solve=='1')
         % iteratively solve all systems simultaneously
            wrk=V(1:n,i);
            tol_solve=1e-4;         % tolerance (change if desired) 
            [v,its_cg]=pcg_prod(K(2:nterm),repmat(wrk,1,nterm-1),100,tol_solve,M1,M1',snew(ir));
         else
         % direct
            wrk=M1*V(1:n,i);
            for kk=2:nterm,
               v(1:n,kk)=(K{kk}+snew(ir)*K{1})\wrk;
            end;
            v=M1'*v;
         end
    else   % polynomial space
        v=multi_matvec(K,M1,M1',V(:,i));
    end
    v = v - V*(V'*v); v = v - V*(V'*v);
    
    % Deselect new basis vectors
    [uu,ss,vv]=svd(v,0);
    ss=diag(ss);
    iv=size(V,2);
    iv_vec(i)=iv;
    if ss(1,1)>1e-12
        addv=1;
        l=cumsum(ss)/sum(ss); il=find(l>=tol_drop,1);
        vnew=uu(:,1:il);
        Mvnew=M1'\vnew;
        V(1:n,iv+1:iv+il)=vnew;  
    else
       addv=0;
    end

    iv=size(V,2);
    iw=size(W,2);
    vector_length=iv*iw;
    
    % Expand projected matrices
    Km{1}=speye(iv);
    ivnew=size(vnew,2);

    for ind=2:nterm
        wrk=K{ind}*Mvnew;
        Mwrk=M1\wrk;
        newk=V(1:n,1:(iv-ivnew))'*Mwrk;
        Km{ind}(1:iv-ivnew,iv-ivnew+1:iv)=newk;  
        Km{ind}(iv-ivnew+1:iv,1:iv-ivnew)=newk';
        Km{ind}(iv-ivnew+1:iv,iv-ivnew+1:iv)=Mwrk'*V(:,iv-ivnew+1:iv);
    end

    if addv, rhs1m=[rhs1m; vnew'*Mrhs1];end
    % Periodically compute the approx solution and residual
    if (rem(tot_it,compute_period)==0)
        rhs2m=rhs2;
        tol_inner=nrmres_noprec*1e-3;
        y0=zeros(iv,iw);
        if (nofirst), 
            y0(1:size(Y,1),1:size(Y,2))=Y;
        end
        [Y,iteraY]=cgkron_m(Km,Gm,rhs1m*rhs2m',y0,iv*iw,tol_inner,iv,iw); % inner solver: cg
        tot_inner=[tot_inner,iteraY];
        nofirst=1;
        % Compute relative variation in the solution
        Yer=Y; [n0,m0]=size(Y0); Yer(1:n0,1:m0)=Yer(1:n0,1:m0)-Y0;
        ErrY=norm(Yer,'fro'); %/normR0;
        nrmres_noprec=ErrY/norm(Y,'fro'); nrmres=nrmres_noprec;
        Y0=Y;
        % Print progress to screen
        fprintf('\n        %2d        %3d       %d      %8.6e         %2d \n', [i,iv,iw,nrmres,iteraY])
    end
    % Save the iteration infos
    stats = savestats();
    info(i+1) = stats;
   error_vec(i+1)=nrmres_noprec;
end

fprintf('\n Total iterations: %d \n\n', i)

%------------------------------------------------
%
% Optional - Estimate Rank + Compute Factors
% 
%-----------------------------------------------

[uu,ss,vv]=svd(Y,0);
ns=sum(diag(ss)/ss(1,1)>tol/iv);
% Determine the solution factors
if (size(V,2)<size(W,1))
   X1=M1'\(V*uu(:,1:ns)); X2=(ss(1:ns,1:ns)*vv(:,1:ns)'*W')';
else
   X1=M1'\(V*(uu(:,1:ns)*ss(1:ns,1:ns))); X2=W*vv(:,1:ns);
end

dimV=size(V,2); final_err=nrmres; avg_inner=mean(tot_inner);
info = info(1:i+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%% ADDED %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Routine in charge of collecting the current iteration stats.
function stats = savestats()
    stats.iter = i;
    % Time
    if i == 0
        stats.time = toc(timetic);
    else
        stats.time = info(i).time + toc(timetic);
    end

    % Ranks and relative residual
    if i == 0
        stats.rank = 0;
        stats.res_norm = rhs_norm;
        stats.res_norm_rel = 1;
    else
        [X1, X2] = get_sol(); X.U = X1; X.V = X2;
        stats.rank = size(X1, 2);
        stats.res_norm = sylv_stable_residual(K, G, C, X);
        stats.res_norm_rel = stats.res_norm / rhs_norm;
    end

    % Restart time
    timetic = tic;
end
function [X1, X2] = get_sol()
    [uu,ss,vv]=svd(Y,0);
    ns=sum(diag(ss)/ss(1,1)>tol/iv);
    % Determine the solution factors
    if (size(V,2)<size(W,1))
       X1=M1'\(V*uu(:,1:ns)); X2=(ss(1:ns,1:ns)*vv(:,1:ns)'*W')';
    else
       X1=M1'\(V*(uu(:,1:ns)*ss(1:ns,1:ns))); X2=W*vv(:,1:ns);
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,k]=pcg_prod(a,b,maxit,tol,L,U,snew)
%function [x,k]=pcg_prod(a,b,maxit,tol,L,U,snew)
% zero initial guess

n=size(b,1);
r=b;
gamma=sum(r.*r); % gamma=r.'*r;
res0=sqrt(gamma);
res=res0;
l=length(a);
x=zeros(n,l);
wrk=x;
mem=ones(l,1);
k=0;

while (max(res./res0) > tol && k<maxit)
    z=r;  
    k=k+1;
    gamma=sum(r.*z);     %gamma=r.'*z;
    if (k==1), p=z;
    else beta=gamma./gamma0;p=z+p*diag(beta);
    end
    ap=U\p;
    for j=1:l
        wrk(1:n,j)=a{j}*ap(1:n,j);
    end
    ap=L\wrk+p*snew;
    delta=sum(p.*ap);
    alfa = gamma./delta;
    x = x + p*diag(alfa);
    r = r - ap*diag(alfa);
    gamma0=gamma;
    res=sqrt(sum(r.*r));
    %  mem=[mem,(res./res0)'];
    %disp([k,res/res0])
end

%max(res./res0)
%fprintf('CG its %d  rel.res %d\n',k,res/res0)


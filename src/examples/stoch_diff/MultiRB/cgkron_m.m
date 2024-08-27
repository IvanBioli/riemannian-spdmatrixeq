function [X,itera]=cgkron_m(K,G,B,X0,maxit,tol,nn,mm)
%function X=cgkron_m(K,G,b,x0,maxit,tol)
% conjugate gradients for hermitian matrices

[n,m]=size(B);
X=X0;

Ax0=stoch_matvec_m(X0,G,K); %AP = a*P; 
R = B-Ax0; 
norma=norm(R,'fro');
res_init=norma;
norma=1;
errore1=1;
errore=errore1;
S = norm(R,'fro')^2;  
P=R; 
itera=0;
beta=0;
tt=norm(P,'fro')^2;

while ( norma>tol & itera<maxit )

  itera=itera+1;
  R1R = S; 
  AP=stoch_matvec_m(P,G,K); %AP = a*P; 
  alfa=R1R/(P(:)'*AP(:)); % alfa=R1R/(P.'*AP); 

  X = X + P * alfa;
  R = R - AP * alfa;
  S = norm(R,'fro')^2; % S=R.'*R; 

  beta=S/R1R;
  P=R+P*beta;
  norma=norm(R(:))/res_init;  
  errore=[errore, norma];

end;


return


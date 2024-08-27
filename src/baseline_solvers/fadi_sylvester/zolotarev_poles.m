function [p, q] = zolotarev_poles(J, a, b, c, d)
% zolotarev_poles - Computes the ADI shifts for AX+XB=F when the 
% eigenvalues of A are in [a,b] and the ones of B in [c,d]

if ~exist('c', 'var')
    c = a;
    d = b;
end

t = c; c = -d; d = -t;

gam = (c-a)*(d-b)/(c-b)/(d-a);                 % Cross-ratio of a,b,c,d
% Calculate Mobius transform T:{-alp,-1,1,alp}->{a,b,c,d} for some alp:
alp = -1 + 2*gam + 2*sqrt(gam^2-gam);          % Mobius exists with this t
A = det([-a*alp a 1; -b b 1 ; c c 1]);         % Determinant formulae for Mobius
B = det([-a*alp -alp a; -b -1 b ; c 1 c]);
C = det([-alp a 1; -1 b 1 ; 1 c 1]);
D = det([-a*alp -alp 1; -b -1 1; c 1 1]);
T = @(z) (A*z+B)./(C*z+D);                     % Mobius transfom

if ( alp < 1e7 )
    K = ellipke( 1-1/alp^2 );                  % ADI shifts for [-1,-1/t]&[1/t,1]
    [~, ~, dn] = ellipj((1/2:J-1/2)*K/J,1-1/alp^2);
else                                           % Prevent underflow when alp large
    K = (2*log(2)+log(alp)) + (-1+2*log(2)+log(alp))/alp^2/4;
    m1 = 1/alp^2; 
    u = (1/2:J-1/2)*K/J; 
    dn = sech(u) + .25*m1*(sinh(u).*cosh(u)+u).*tanh(u).*sech(u); 
end
p = T( -alp*dn ); q = T( alp*dn );             % ADI shifts for [a,b]&[c,d]
end

function omega = get_nodes2(m,M,ell)
%function omega = get_nodes2(m,M,ell)
%
% Optimal parameters for the given interval and ell nodes
% reference. p. 861 in Alternating Direction Implicit iteration for systems
% with complex spectra}}, Nancy~S. Ellner and Eugene~L. Wachspress},
% SIAM J. Numer. Anal, 23 (3), (1991), pp.{859--870}.
% or  "The ADI Model Problem", E. Wachspress (book) 
%
% For the code, see, e.g., J.Sabino, Solution of Large-Scale Lyapunov 
% Equations via  the Block Modified Smith Method,
% PhD Thesis, Rice University, 2006. p.43
% A good approximation is given by logarithmically distributed
% points in [a,b]:    omega = logspace(log10(a),log10(b),n)’;
%

kk=1-(m/M)^2;
[K,E]=ellipke(kk);
t = (2*(1:ell)-1)*K/(2*ell);
[sn,cn,dn]=ellipj(t,kk);
omega =sort(dn*M)';

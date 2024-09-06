% drawn random beta binomial numbers

%drawn n samples with a
%probability distribution with a mean of p, and shape a (corresponds to
%parameter alpha in beta distribution)

function R= betabinornd(n,p,rho)%,M,N)
if n == 0
    R = 0;
else
    a = p.*(1/rho - 1);
    b=a.*(1-p)./p;


% if nargin<4
    
    R=binornd(n,betarnd(a,b));
    
% else
    
%     R=binornd(n,betarnd(a,b,M,N));

end
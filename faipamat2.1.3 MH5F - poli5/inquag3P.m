%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function r=inquag3P(g1,g2,g3,gtarg,t2,t3)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   FAIPAMAT 2.1 - 23/12/1999
%
%   QUADRATIC interpolation of the constraints to find an interior point
%
%   BASED ON 3 POINTS without derivatives
%
%   t1=0; t2 and t3


c=g1;
deps = 1.e-10; 
[coef] = [t2^2 t2; t3^2 t3 ]\[g2-g1;g3-g1];
a=coef(1); b=coef(2);
pcoef = [a; b; c - gtarg];
rv = roots(pcoef);
nroots = length(rv);
if nroots == 0     
        r = inf;
else
        for k = 1:nroots
           if imag(rv(k)) ~= 0, rv(k)=inf; end
           if rv(k) < deps, rv(k)=inf; end 
        end
        r = min(rv);
end



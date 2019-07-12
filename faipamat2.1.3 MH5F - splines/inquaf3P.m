%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function r=inquaf3P(f1,f2,f3,t2,t3)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   FAIPAMAT 2.1 - 23/12/1999
%
%   QUADRATIC interpolation of the function
%
%   BASED ON 3 POINTS without derivatives:
%
%   t1=0; t2 and t3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c=f1;
deps = 1.e-10; %
[coef] = [t2^2 t2; t3^2 t3 ]\[f2-f1;f3-f1];
a=coef(1); b=coef(2);
r = -b/(2*a);
if r <= 0 , r= inf; end  % 23/12/99
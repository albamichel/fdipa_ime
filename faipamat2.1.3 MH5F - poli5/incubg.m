%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function r=incubg(fnew,fold,ftarg,graddnew,graddold,stepsize)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   FAIPAMAT 2.0 - 27/09/1999
%
%
%   CUBIC interpolation of the constraints to find an interior point
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deps = 1.e-10; 
u = (fnew - fold)/stepsize;
v = graddnew + graddold - 3*u; 
a = 3*(graddnew + graddold - 2*u)/(stepsize^2);
%if a < deps, a = 0; end 
b = -(graddold + v)/stepsize;
%if a ==0 & b < deps, b=0; end 
pcoef = [a/3; b; graddold; fold-ftarg];
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

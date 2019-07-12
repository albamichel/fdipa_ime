%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function r=inquadf2P(fnew,fold,gradold,t)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   FAIPAMAT 2.0 - 27/09/1999
%
%   Quadratic interpolation based on 2 points and their gradients to 
%   estimate the minimum.
%
%   Uses the function and the derivative at t=0 and the derivative
%   at t=stepsize
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deps = 1.E-10;
b = gradold;
c = fold;
a = (fnew - t*b - c)/t^2;
r = -b/(2*a);

if r < deps
   r = inf;
end 
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function r=incubf(fnew,fold,graddnew,graddold,stepsize)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   FAIPAMAT 2.0 - 27/09/1999
%
%
%   This function uses cubic interpolation and the values of two 
%   points and their gradients in order to estimate the minimum of
%   a function along a line.
%   Formulation according C. Lemarechal

deps = 1.E-10; 
u = (fnew - fold)/stepsize;
v = graddnew + graddold - 3*u; 
a = 3*(graddnew + graddold - 2*u)/(stepsize^2);
b = -(graddold + v)/stepsize;
delta = (v*v - graddnew*graddold)/(stepsize^2); 
        % delta = b*b - a* graddold
%if ~(delta > 0) ...  % Changed 27/05/99
if (delta < 0) ...  % The cubic polinome has not a minimum with t > 0
        r = inf;
else
   if b > 0
      r = -graddold/(b + sqrt(delta));
   else
      if abs(a) < deps % The interpolation polinome is almost
                       % quadratic, with negative curvature.
            r = inf;
      else
            r = (-b + sqrt(delta))/a;
      end
   end
end

if r < deps

%   disp('WARNING, negative root in incubf- takes r = inf');
   
%  TAKES r = inf when a very little or negative step
%        was obtained

    r = inf;
end 
 

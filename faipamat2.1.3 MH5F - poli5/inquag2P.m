%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function r=inquag2P(fnew,fold,ftarg,gradold,t)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   FAIPAMAT 2.0 - 27/09/1999
%
%   QUADRATIC interpolation of the constraints to find an interior point
%
%   BASED ON 2 POINTS  and the derivative at the left
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deps = 1.e-6;
 
b = gradold;
c  = fold; 
a = (fnew - b*t - c)/(t^2);

if abs(a) < deps              
   a=0;
end  

if a==0
 
% Linear interpolation 

   a = (fnew - fold)/t;
   r=(ftarg -fold)/a;
        
else 
   delta = b^2-4*a*(c-ftarg);
   if delta <0
      r=inf;
   else
      r1=(-b-delta^.5)/(2*a);
      r2=(-b+delta^.5)/(2*a);
      r=minpos([r1 r2 inf ]);
   end   
end      

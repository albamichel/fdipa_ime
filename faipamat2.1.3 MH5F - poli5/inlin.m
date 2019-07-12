%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function r=inlin(fnew,fold,ftarg,t)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   FAIPAMAT 2.0 - 27/09/1999
%
%
%   LINEAR interpolation of the constraints to find an interior point.
%
%   BASED ON the value of the constraints at 2 POINTS.
% 

a = (fnew - fold)/t;

r=(ftarg -fold)/a; 

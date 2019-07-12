%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function[]=dispstop(istop)
%
%
%   FAIPAMAT 2.0 - 24/09/1999
% 
%   PRINTS STOPPING CRITERIUM, WHEN IT STOPS

if istop~=0
disp('                                                             ')
disp('=============================================================')
disp('                                                             ')
disp('STOPPING CRITERION:')
disp('                                                             ')
end

if istop==1
disp('Gradient of the Lagrangian less than tolerance')
end

if istop==2
disp('Descent direction "d0" less than tolerance')
end

if istop==3
disp('Reduction of the penalty function less than tolerance')
end

if istop==4
disp('Maximum number of main iterations')
end

if istop==5
disp('Maximum number of iterations in line search')
end

if istop==6
disp('The objective function is lower than the expected given value')
end

if istop==10 | istop==11
disp('A Feasible point was obtained')
end
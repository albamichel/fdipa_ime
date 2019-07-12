%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function[istop]=stoptest(istop,modd0,mglag0,f,g_all,oldpen,pen,erreq,iter,...
                         iterlin,data,idata)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   FAIPAMAT 2.0 - 27/09/1999
%
%   Verifies stopping criteria
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  

ifpoint  = idata(15);
                                      
if mglag0 < data(1),
   istop=1; % Gradient of the Lagrangian less than tolerance
elseif modd0 < data(2),
   istop=2; % Descent direction "d0" less than tolerance
elseif (oldpen - pen) < data(3) & oldpen > pen 
   istop=3; % Reduction of the penalty function less than tolerance
elseif f < data(5) &  ifpoint~=1; % The function is lower than the 
                                 % expected given value
   istop=6;
elseif ifpoint~=0 & all(g_all+f < - eps)
% elseif ifpoint~=0 & all(g_all+f < 0) Modificado em 24/09/2001
   
   istop=10;       
end
if ifpoint~=0 & (istop==1 | istop==2 | istop==3)
   disp('WARNING:  A FEASIBLE POINT WAS NOT OBTAINED ')
end   
if ifpoint==2 & istop==10, istop=11; end

if erreq > data(4),istop=0;end     % Tests error of the equality const.

if iter < idata(8) & ifpoint==0, istop=0; end   % Less iterations that the minimum.

if iter == idata(7), istop=4; end  % Maximum number of main iterations.

if iterlin >= idata(9),            % Maximum number of iterations in
   istop=5;                        % line search.                       
end                                  

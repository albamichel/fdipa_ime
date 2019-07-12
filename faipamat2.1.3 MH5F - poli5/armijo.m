%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [t,istalin,cod]=armijo(t,f,g_all,oldf,oldgdf,neqlin,ncstr,data);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   FAIPAMAT 2.0 - 23/09/1999
% 
%   ARMIJO'S LINE SEARCH 
               
istalin=0;

% PARAMETERS FOR STOPPING CRITERIA IN LINE SERACH: 
              
eta1 = data(16);
nu   = data(18);
deps = 1.e-7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


if f <= oldf+t*eta1*oldgdf & all(g_all(neqlin+1:ncstr)<=0) 
   
   cod = 2;      % ARMIJO'S CRITERIUM IS VERIFIED
   istalin = 1;
   
else % (f <= oldf+t*eta1*oldgdf & all(g_all(neqlin+1:ncstr)<=0))
   
   cod = 1; % Defines a new step
   t = nu*t; 
   
end

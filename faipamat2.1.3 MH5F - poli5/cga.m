

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [var] = cga(A,var,b,errcg1,errcg2);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   FAIPAMAT 2.0 - 23/09/1999
% 
%   Solves a linear symmetric system of equations     
%   with the Congugate Gradient Algorithm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

d = b - (A*var);
g = -d;
k = 0;
normb=norm(b);
if normb==0, normb=1; end
while max(abs(g/normb))>errcg2 & max(abs(g)) > errcg1,
   k = k+1;
   alfa = -(g'*d)/(d'*A*d);
   var = var + alfa*d;                        
   g = (A*var) - b;                        
   beta = (g'*A*d)/(d'*A*d);             
   d = -g + (beta*d);                     
end;

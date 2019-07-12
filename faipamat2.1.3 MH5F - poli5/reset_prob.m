%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function[x,g,lvlb,lvub,nvar,ncstr,neq,neqlin] = ...
...
reset_prob(x,g,lvlb,lvub,nvar,ncstr,neq,neqlin,neq_op,neqlin_op)
%                                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Reset to the original problem, after obtaining a feasible point
%
%   FAIPAMAT 2.1 - 23/12/1999
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


z = x(nvar);
nvar  = nvar - 1;
neq = neq_op;
ncstr = ncstr + neq;
neqlin = neqlin_op;
x = x(1:nvar);
g = g + z;
lvlb = lvlb(1:nvar); 
lvub = lvub(1:nvar);
                                
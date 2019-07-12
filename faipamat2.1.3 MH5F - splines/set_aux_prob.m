%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function[x,g,lvlb,lvub,nvar,ncstr,ntcstr,neq,neqlin,B,DDELTA,GGAMA,...
         d2,glag2,omegaE,c] = ...
...
set_aux_prob(x,g,lvlb,lvub,nvar,ncstr,ntcstr,neq,neqlin,B,DDELTA,GGAMA,...
             d2,glag2,omegaE,c,idata)
%                                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Defines the auxiliary problem:
%
%          Find (x,z) that minimizes z
%
%                             s.t. g(x) <= z
%
% that gives a feasible point when z <0.
%
%
%   FAIPAMAT 2.1 - 23/12/1999
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


iupdateb = idata(4);
nvar  = nvar + 1;
ncstr = ncstr - neq;
ntcstr=ntcstr-neq;
neq = 0;
neqlin = 0;
z = max(1.2* max(g),.00001); %>>>>>>>>>>>>> 30/8/200
x = [x;z];
g = g - z;
lvlb = [lvlb(:); 0];
lvub = [lvub(:); 0];
if iupdateb == 1 | iupdateb == 0
   B = [B zeros(nvar-1,1);zeros(1,nvar-1) 1];
end

DDELTA=zeros(nvar,0); GGAMA=zeros(nvar,0); 
d2=zeros(nvar,1);glag2=zeros(nvar,1);
omegaE=zeros(0,1);
c = zeros(neq,1);
                                
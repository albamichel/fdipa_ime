%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function[lambda0]=lam_box(nvar,nbind,d0,lambda0,lambda,glow,gup,...
                               lvlb,lvub);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   FAIPAMAT 2.0 - 27/09/1999
%                               
%   Calculus of the lagrange multipliers corresponding to the box 
%   constraints: 
%
%   lambda0=-G_1^(-1)*diag(lambda)*gradg'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
      
j=0;
for i=1:nvar
   if lvlb(i)==1 
      j=j+1;
      lambda0 = [ lambda0; 1/glow(j)*lambda(nbind+j)*d0(i)];
   end
end
lenvlb=j;
j=0;
for i=1:nvar
   if lvub(i)==1 
      j=j+1;
      lambda0=[ lambda0; -1/gup(j)*lambda(nbind+lenvlb+j)*d0(i)];
   end
end

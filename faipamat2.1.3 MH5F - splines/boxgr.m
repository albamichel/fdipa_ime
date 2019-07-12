%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[gradbox]=boxgr(nvar,lenvlb,lenvub,lvlb,lvub)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   COMPUTES THE MATRICES OF BOX CONSTRAINTS DERIVATIVES
%
%
%   FAIPAMAT 15/11/2000
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      gradbox=zeros(nvar,lenvlb+lenvub);
      j=0;
      jj=lenvlb;
      for i=1:nvar
         if lvlb(i)==1 
            j=j+1;          
            gradbox(i,j)=-1; 
         end
         if lvub(i)==1
            jj=jj+1;
            gradbox(i,jj)=1; 
         end
      end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[gvlb,gvub]=boxderiv(nvar,lenvlb,lenvub,lvlb,lvub)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   COMPUTES THE MATRICES OF BOX CONSTRAINTS DERIVATIVES
%
%
%   FAIPAMAT 2.0 - 24/09/1999
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      gvlb=zeros(nvar,lenvlb);
      gvub=zeros(nvar,lenvub);
      j=0;
      jj=0;
      for i=1:nvar
         if lvlb(i)==1 
            j=j+1;          
            gvlb(i,j)=-1; 
         end
         if lvub(i)==1
            jj=jj+1;
            gvub(i,jj)=1; 
         end
      end
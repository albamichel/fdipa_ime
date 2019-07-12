%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[gboxM]=boxgrxM(nvar,lenvlb,lenvub,lvlb,lvub,M,colM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   COMPUTES THE MATRICES OF BOX CONSTRAINTS DERIVATIVES
%
%
%   FAIPAMAT 15/11/2000
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     M=M(lenvlb+lenvub,colM)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      gboxM=zeros(nvar,colM);

      for k=1:colM
         j=0;
         jj=lenvlb;
         for i=1:nvar
            if lvlb(i)==1 
               j=j+1;          
               gboxM(i,k) = gboxM(i,k) -  M(j,k); 
            end
            if lvub(i)==1
               jj=jj+1;
               gboxM(i,k)= gboxM(i,k) + M(jj,k); 
            end
         end
      end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function MHM = MxHxM(M,DELTA,GAMA,RMAT,DMAT,GAMAtGAMA,npairs);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   FAIPAMAT 2.0 - 23/09/1999
% 
%
%   Computes the product Mt*H*M where H is a limited memory 
%
%   quasi-Newton Matrix


linM = size(M,1);
colM = size(M,2);
MHM = zeros(colM,colM);
for j=1:colM
   MHM(:,j) = M'*Hxv(M(:,j),DELTA,GAMA,RMAT,DMAT,GAMAtGAMA,npairs);
end   


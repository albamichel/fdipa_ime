%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Hv = Hxv(v,DELTA,GAMA,RMAT,DMAT,GAMAtGAMA,npairs)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   FAIPAMAT 2.0 - 28/09/1999
% 
%*******************************************************************************%
%     Computes the product of a Quasi-Newton Matrix obtained by limited         %
%     memory by a vector: Hv=h*v                                                %
%                                                                               %
%        v                => vetor to be multiplied at the left by H            %
%        DELTA, GAMA, RMAT, DMAT, GAMAtGAMA => auxiliar matrices                %
%        npairs                => number of stored pairs                        %
%                                                                               %
%     Ref.: Byrd, R. H., Nocedal, J. and Schnabel, R. B., "Representation       %
%     of Quasi-Newton Matrices and Their Use in Limited Memory Methods",        %
%     Technical Report CU-CS-612-92, 1992, University of Colorado at Boulder    %
%*******************************************************************************%

if npairs==0
   Hv = v;
else
%   Normaliztion
%
%   gama = 1 /(RMAT(npairs,npairs)*GAMAtGAMA(npairs,npairs));
%
   gama=1;
   A = RMAT * (DELTA' * v);
   p = [RMAT'*((spdiags(DMAT',0,npairs,npairs)+gama*GAMAtGAMA)*A-gama*GAMA'*v);-A];
   Hv = gama*v+[DELTA,gama*GAMA]*p;
end
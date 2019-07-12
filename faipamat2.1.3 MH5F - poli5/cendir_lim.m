%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function[d2,glag2] = cendir_lim(DELTA,GAMA,RMAT,DMAT,GAMAtGAMA,npairs,gt,...
                                gradg,lvlb,lvub,lambda,omegaI2,omegaE2,...
                                neq,nvar,nbind,lenvlb,lenvub,isolver);
%                                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   FAIPAMAT 2.0 - 24/09/1999
%          
%
% COMPUTES THE CENTERING DIRECTION "d2", WITH LIMITED MEMORY QUASI-NEWTON         

  ntbind=neq+nbind+lenvlb+lenvub;
  

% Derivatives of the box constraints  
  
  
  [gvlb,gvub]=boxderiv(nvar,lenvlb,lenvub,lvlb,lvub);


% Computing d2


if isolver == 2

  % By solving the nonsymmetric system in (lambda,mu)
    
   mat = MxHxM([gradg gvlb gvub],DELTA,GAMA,RMAT,DMAT,GAMAtGAMA,npairs);
   mat = diag([ones(neq,1);lambda])*mat;
   mat(neq+1:ntbind,neq+1:ntbind) =  mat(neq+1:ntbind,neq+1:ntbind) ...
                                    - diag(gt(neq+1:ntbind));
    
   mataux3 = mat\[ omegaE2;diag(lambda(1:nbind))*omegaI2;zeros(lenvlb+lenvub,1)];               
          
   lambda2 = mataux3(neq+1:ntbind,1);                                                    
   mu2 = mataux3(1:neq,1);
 
   clear mataux3
   glag2 = ([gradg(:,neq+1:neq+nbind) gvlb gvub]*lambda2 ...
             + gradg(:,1:neq)*mu2);

   d2 = - Hxv(glag2,DELTA,GAMA,RMAT,DMAT,GAMAtGAMA,npairs);      
             
elseif isolver == 3

 % By solving the symmetric system in (lambda,mu)
 
 
   mat = MxHxM([gradg gvlb gvub],DELTA,GAMA,RMAT,DMAT,GAMAtGAMA,npairs);              
   mat(neq+1:ntbind,neq+1:ntbind) =  mat(neq+1:ntbind,neq+1:ntbind) ...
                                    - (diag(gt(neq+1:ntbind)./lambda));
                                    
   mataux3 = mat\[ omegaE2;omegaI2;zeros(lenvlb+lenvub,1)];
 
   lambda2 = mataux3(neq+1:ntbind,1);                                                    
   mu2 = mataux3(1:neq,1); 
   clear mataux3
   glag2 = ([gradg(:,neq+1:neq+nbind) gvlb gvub]*lambda2 ...
             + gradg(:,1:neq)*mu2);             
   d2 = - Hxv(glag2,DELTA,GAMA,RMAT,DMAT,GAMAtGAMA,npairs); 
             
end

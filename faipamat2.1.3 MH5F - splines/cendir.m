%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function[d2,glag2] = cendir(B,gt,gradg,lvlb,lvub,lambda,omegaI2,...
                      omegaE2,neq,nvar,nbind,lenvlb,lenvub,...
                      isolver,errcg1,errcg2);
%                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
%
%   FAIPAMAT 2.0 - 23/09/1999
% 
%         
%   COMPUTES THE CENTERING DIRECTION "d2"         

  ntbind=neq+nbind+lenvlb+lenvub;
  

% Derivatives of the box constraints  
  
if isolver~=1  
  [gvlb,gvub]=boxderiv(nvar,lenvlb,lenvub,lvlb,lvub);
end
% Computing d2

if isolver == 1

% By solving the system in (d,lambda,mu):
   mataux1=zeros(nvar+ntbind,nvar+ntbind);
   mataux1(1:nvar,nvar+neq+nbind+1:nvar+ntbind)= ...
      boxgr(nvar,lenvlb,lenvub,lvlb,lvub);
   
   mataux1 = [B [gradg mataux1(1:nvar,nvar+neq+nbind+1:nvar+ntbind)]; ...  
              gradg(:,1:neq)' zeros(neq,ntbind); ...
              diag(lambda)*[gradg(:,neq+1:neq+nbind) ...
                            mataux1(1:nvar,nvar+neq+nbind+1:nvar+ntbind)]'... 
              zeros(ntbind-neq,neq) diag(gt(neq+1:ntbind))];

  


 
   mataux1 = mataux1 \[ zeros(nvar,1); -omegaE2; ...
             -diag(lambda(1:nbind))*omegaI2; ...
             -zeros(lenvlb+lenvub,1)];
             
   lambda2 = mataux1(nvar+neq+1:nvar+ntbind,1);
   mu2 = mataux1(nvar+1:nvar+neq,1);
   d2 = mataux1(1:nvar,1);
   glag2 = gradg(:,neq+1:neq+nbind)*lambda2(1:nbind,1) ...
             + boxgrxM(nvar,lenvlb,lenvub,lvlb,lvub,...
                      lambda2(nbind+1:nbind+lenvlb+lenvub),1) ...
              + gradg(:,1:neq)*mu2;

   clear mataux1

elseif isolver == 2

  % By solving the nonsymmetric system in (lambda,mu)
                 
   mat = diag([ones(neq,1);lambda])*[gradg gvlb gvub]'*B*[gradg gvlb gvub];
   mat(neq+1:ntbind,neq+1:ntbind) =  mat(neq+1:ntbind,neq+1:ntbind) ...
                                    - diag(gt(neq+1:ntbind));
                  
   mataux3 = mat\[ omegaE2;diag(lambda(1:nbind))*omegaI2;zeros(lenvlb+lenvub,1)];               
          
   lambda2 = mataux3(neq+1:ntbind,1);                                                    
   mu2 = mataux3(1:neq,1);
 
   clear mataux3
   glag2 = ([gradg(:,neq+1:neq+nbind) gvlb gvub]*lambda2 ...
             + gradg(:,1:neq)*mu2);

   d2 = - B*glag2;      
             
elseif isolver == 3

 % By solving the symmetric system in (lambda,mu)
             
   mat = [gradg gvlb gvub]'*B*[gradg gvlb gvub];
   mat(neq+1:ntbind,neq+1:ntbind) =  mat(neq+1:ntbind,neq+1:ntbind) ...
                                    - (diag(gt(neq+1:ntbind)./lambda));
                                    
   mataux3 = mat\[ omegaE2;omegaI2;zeros(lenvlb+lenvub,1)];
 
   lambda2 = mataux3(neq+1:ntbind,1);                                                    
   mu2 = mataux3(1:neq,1);
 
   clear mataux3
   glag2 = ([gradg(:,neq+1:neq+nbind) gvlb gvub]*lambda2 ...
             + gradg(:,1:neq)*mu2);

   d2 = - B*glag2;
           
elseif isolver == 4

  % By eliminating lambda and solving a system in d0 and mu 
  
    d2 = [(B - [gradg(:,neq+1:neq+nbind) gvlb gvub]*...                                               
          diag(lambda./gt(neq+1:ntbind))*...
          [gradg(:,neq+1:neq+nbind) gvlb gvub]') ...
          gradg(:,1:neq); gradg(:,1:neq)' zeros(neq,neq)]...
         \[gradg(:,neq+1:neq+nbind)*...
          diag(lambda(1:nbind)./gt(neq+1:neq+nbind))*omegaI2;...
           -omegaE2];
                         
    d2 = d2(1:nvar,1);
    glag2=[];
    
elseif isolver == 5

 % By solving the symmetric system in (lambda,mu)
             
   mat = [gradg gvlb gvub]'*B*[gradg gvlb gvub];
   mat(neq+1:ntbind,neq+1:ntbind) =  mat(neq+1:ntbind,neq+1:ntbind) ...
                                    - (diag(gt(neq+1:ntbind)./lambda));
                                    
   b1 = [ omegaE2;omegaI2;zeros(lenvlb+lenvub,1)];
   init = zeros(length(mat),1);
   mataux3(:,1) = cga(mat,init,b1,errcg1,errcg2); 
   lambda2 = mataux3(neq+1:ntbind,1);                                                    
   mu2 = mataux3(1:neq,1);
 
   clear mataux3
   glag2 = ([gradg(:,neq+1:neq+nbind) gvlb gvub]*lambda2 ...
             + gradg(:,1:neq)*mu2);

   d2 = - B*glag2;
             
elseif isolver == 6
                         
   % By solving the symmetric system in (lambda,mu) with the
   % CONJUGATE GRADIENT ALGORITHM

   mat = [gradg gvlb gvub]'*B*[gradg gvlb gvub];
   mat(neq+1:ntbind,neq+1:ntbind) =  mat(neq+1:ntbind,neq+1:ntbind) ...
                                    - (diag(gt(neq+1:ntbind)./lambda));
   b = [ omegaE2;omegaI2;zeros(lenvlb+lenvub,1)];
   init = zeros(length(mat),1);
   res = cga(mat,init,b,errcg1,errcg2);                    
   lambda2 = res(neq+1:ntbind,1);                                                    
   mu2 = res(1:neq,1);
 
   clear res
   glag2 = ([gradg(:,neq+1:neq+nbind) gvlb gvub]*lambda2 ...
             + gradg(:,1:neq)*mu2);

   d2 = - B*glag2;
       
end

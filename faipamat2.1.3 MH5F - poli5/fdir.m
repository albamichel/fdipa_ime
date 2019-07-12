%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function[d,d0,d02,lambda0,mu0,c,pen,gpen,glag,mglag0] = fdir(B,f,gt,...
         gradf,gradg,lvlb,lvub,lambda,omegaI,omegaE,c,neq,nvar,...
         nbind,lenvlb,lenvub,isolver,errcg1,errcg2,data);
%         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FAIPAMAT 2.0
%
% COMPUTES THE FEASIBLE DESCENT DIRECTION " d = d0 + rho*d1 "


fi=data(11);
alfa=data(12);
glag0=[];
glag1=[];  
ntbind=neq+nbind+lenvlb+lenvub;
  
% Derivatives of the box constraints

if isolver ~=1
   [gvlb,gvub]=boxderiv(nvar,lenvlb,lenvub,lvlb,lvub);
end
% Computing d0 and d1:
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
        
   mataux1 = mataux1 \[[-gradf; -gt(1:neq) ; zeros(ntbind-neq,1)] ...
                       [ zeros(nvar,1);-omegaE;-diag(lambda)*omegaI]];       
                 
   lambda0 = mataux1(nvar+neq+1:nvar+ntbind,1);
   mu0 = mataux1(nvar+1:nvar+neq,1); 
   lambda1 = mataux1(nvar+neq+1:nvar+ntbind,2);
   mu1 = mataux1(nvar+1:nvar+neq,2); 
   d0 = mataux1(1:nvar,1);
   d1 = mataux1(1:nvar,2);
  
   glag0 =  (gradf + gradg(:,neq+1:neq+nbind)*lambda0(1:nbind,1) ...
             + boxgrxM(nvar,lenvlb,lenvub,lvlb,lvub,...
                      lambda0(nbind+1:nbind+lenvlb+lenvub),1)...
             + gradg(:,1:neq)*mu0);     
   glag1 = gradg(:,neq+1:neq+nbind)*lambda1(1:nbind,1) ...
             + boxgrxM(nvar,lenvlb,lenvub,lvlb,lvub,...
                      lambda1(nbind+1:nbind+lenvlb+lenvub),1) ...
              + gradg(:,1:neq)*mu1;
  
   clear mataux1
   
elseif isolver == 2   
   % By solving the nonsymmetric system in (lambda,mu)
        
   mat = diag([ones(neq,1);lambda])*[gradg gvlb gvub]'*B*[gradg gvlb gvub];
   mat(neq+1:ntbind,neq+1:ntbind) =  mat(neq+1:ntbind,neq+1:ntbind) ...
                                    - diag(gt(neq+1:ntbind));
                  
   mataux3 = mat\[(-diag([ones(neq,1);lambda])*(B*[gradg gvlb gvub])'*gradf + ...
                  [gt(1:neq); zeros(ntbind-neq,1)]) ...
                  [omegaE; diag(lambda)*omegaI]];                                         
   
   lambda0 = mataux3(neq+1:ntbind,1);
   mu0 = mataux3(1:neq,1);           
   lambda1 = mataux3(neq+1:ntbind,2);
   mu1 = mataux3(1:neq,2); 
   
   glag0 =  (gradf + [gradg(:,neq+1:neq+nbind) gvlb gvub]*lambda0 ...
             + gradg(:,1:neq)*mu0);     
   d0 = - B*glag0;
   glag1 = ([gradg(:,neq+1:neq+nbind) gvlb gvub]*lambda1 ...
              + gradg(:,1:neq)*mu1);       
   d1 = - B*glag1;
   clear mataux3
   clear mat                   
   
elseif isolver == 3                         
   % By solving the symmetric system in (lambda,mu)
   
   mat = [gradg gvlb gvub]'*B*[gradg gvlb gvub];
   mat(neq+1:ntbind,neq+1:ntbind) =  mat(neq+1:ntbind,neq+1:ntbind) ...
                                    - (diag(gt(neq+1:ntbind)./lambda));
                  
   mataux3 = mat\[(-(B*[gradg gvlb gvub])'*gradf + ...
                  [gt(1:neq); zeros(ntbind-neq,1)]) ...
                  [omegaE; omegaI]];                  
              
   lambda0 = mataux3(neq+1:ntbind,1);
   mu0 = mataux3(1:neq,1);           
   lambda1 = mataux3(neq+1:ntbind,2);
   mu1 = mataux3(1:neq,2); 
   
   glag0 =  (gradf + [gradg(:,neq+1:neq+nbind) gvlb gvub]*lambda0 ...
             + gradg(:,1:neq)*mu0);     
   d0 = - B*glag0;
   glag1 = ([gradg(:,neq+1:neq+nbind) gvlb gvub]*lambda1 ...
              + gradg(:,1:neq)*mu1);       
   d1 = - B*glag1;
   clear mataux3
   clear mat              

elseif isolver == 4                         
   % By eliminating lambda and solving a system in d0 and mu
    
   mataux3 = [(B - [gradg(:,neq+1:neq+nbind) gvlb gvub]*...                                               
                         diag(lambda./gt(neq+1:ntbind))*...
                         [gradg(:,neq+1:neq+nbind) gvlb gvub]') ...
              gradg(:,1:neq); gradg(:,1:neq)' zeros(neq,neq)]...
             \[[-gradf; -gt(1:neq)] ...
              [[gradg(:,neq+1:neq+nbind) gvlb gvub]*...                                               
              diag(lambda./gt(neq+1:ntbind))*omegaI; -omegaE]];
                         
   d0 = mataux3(1:nvar,1);
   d1 = mataux3(1:nvar,2);
   mu0 = mataux3(nvar+1:nvar+neq,1);                                                    
   mu1 = mataux3(nvar+1:nvar+neq,2);
   lambda0 = - diag(lambda./gt(neq+1:ntbind))*...
               [gradg(:,neq+1:neq+nbind) gvlb gvub]'*d0;
   lambda1 = - diag(lambda./gt(neq+1:ntbind))*...
               [gradg(:,neq+1:neq+nbind) gvlb gvub]'*d1...
                - diag(lambda./gt(neq+1:ntbind))*omegaI;
   clear mataux3
          
   
elseif isolver == 5
                         
   % By solving the symmetric system in (lambda,mu) with the
   % CONJUGATE GRADIENT ALGORITHM

   mat = [gradg gvlb gvub]'*B*[gradg gvlb gvub];
   mat(neq+1:ntbind,neq+1:ntbind) =  mat(neq+1:ntbind,neq+1:ntbind) ...
                                    - (diag(gt(neq+1:ntbind)./lambda));
                  
   b1 =(-(B*[gradg gvlb gvub])'*gradf + ...
       [gt(1:neq); zeros(ntbind-neq,1)]); 
   b2 = [omegaE; omegaI];                  
   init = zeros(length(mat),1);
   res1 = cga(mat,init,b1,errcg1,errcg2); 
   lambda0 = res1(neq+1:ntbind,1); 
   mu0 = res1(1:neq,1);
   glag0 =  (gradf + [gradg(:,neq+1:neq+nbind) gvlb gvub]*lambda0 ...
             + gradg(:,1:neq)*mu0);     
   d0 = - B*glag0;
   d02 = d0'*d0;
%   c = max(c,-1.2*mu0); MODIFICADO EM 4/09/2003 
   c = max(c,abs(1.2*mu0));

   pen = f;
   gpen = gradf;
   if neq > 0
      pen = pen - gt(1:neq)'*c;
      gpen = gpen - gradg(:,1:neq)*c;
   end

   if [lambda0;mu0+c]'*[omegaI;omegaE] > 0,
      ro = min(fi*d02,(alfa-1)*(d0'*gpen)/([lambda0;mu0+c]'*[omegaI;omegaE]));
   else
      ro = fi*d02;
   end
   b1 = -(B*[gradg gvlb gvub])'*gradf + ...
       [zeros(neq,1) + (diag(c./(mu0+c)))*gt(1:neq); zeros(ntbind-neq,1)];
   mat(1:neq,1:neq) = mat(1:neq,1:neq) - diag(gt(1:neq)./(mu0+c));        
   res2 = cga(mat,res1,b1+ro*b2,errcg1,errcg2);  

   lambda = res2(neq+1:ntbind,1);
   mu = res2(1:neq,1); 
   
   glag = (gradf + [gradg(:,neq+1:neq+nbind) gvlb gvub]*lambda ...
             + gradg(:,1:neq)*mu); 
   d = - B * glag;

   clear mat res1 res1 b1 b2
elseif isolver==6

% By solving the symmetric system in (lambda,mu)
   
   mat = [gradg gvlb gvub]'*B*[gradg gvlb gvub];
   mat(neq+1:ntbind,neq+1:ntbind) =  mat(neq+1:ntbind,neq+1:ntbind) ...
                                    - (diag(gt(neq+1:ntbind)./lambda));
   b1 =(-(B*[gradg gvlb gvub])'*gradf + ...
       [gt(1:neq); zeros(ntbind-neq,1)]); 
   b2 = [omegaE; omegaI];                  
   init = zeros(length(mat),1);
   mataux3(:,1) = cga(mat,init,b1,errcg1,errcg2);                                     
   mataux3(:,2) = cga(mat,init,b2,errcg1,errcg2);             
   lambda0 = mataux3(neq+1:ntbind,1);
   mu0 = mataux3(1:neq,1);           
   lambda1 = mataux3(neq+1:ntbind,2);
   mu1 = mataux3(1:neq,2); 
   
   glag0 =  (gradf + [gradg(:,neq+1:neq+nbind) gvlb gvub]*lambda0 ...
             + gradg(:,1:neq)*mu0);     
   d0 = - B*glag0;
   glag1 = ([gradg(:,neq+1:neq+nbind) gvlb gvub]*lambda1 ...
              + gradg(:,1:neq)*mu1);       
   d1 = - B*glag1;
   clear mataux3
   clear mat         
   
end   

if isolver ~= 5,
   % Updating c and the penalty function and its gradient
%   c = max(c,-1.2*mu0); MODIFICADO EM 4/09/2003 
   c = max(c,abs(1.2*mu0));   pen = f;
   gpen = gradf;
   if neq > 0
      pen = pen - gt(1:neq)'*c;
      gpen = gpen - gradg(:,1:neq)*c;
   end

   % Updating ro

   d02 = d0'*d0;
   if d1'*gpen > 0
      ro = min(fi*d02, (alfa-1)*(d0'*gpen)/(d1'*gpen));
   else
      ro = fi*d02;
   end
        
   % The search direction
       
   d = d0 + ro*d1;
   % Gradient of the lagrangian corresponding to 
   % lambda_bar=lamda0+ro*lambda1

   glag = glag0 + ro*glag1;
   
end   

% Lagrangian's gradient norm
%if isolver==3
   mglag0=norm(glag0);                     
%else
%   mglag0=norm(gradf+[gradg gvlb gvub]*[mu0;lambda0]); 
%end

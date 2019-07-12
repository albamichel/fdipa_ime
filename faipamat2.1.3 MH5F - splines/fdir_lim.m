%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function[d,d0,d02,lambda0,mu0,c,pen,gpen,glag,mglag0] = ...
         fdir_lim(DELTA,GAMA,RMAT,DMAT,GAMAtGAMA,npairs,f,gt,...
         gradf,gradg,lvlb,lvub,lambda,omegaI,omegaE,c,neq,nvar,...
         nbind,lenvlb,lenvub,isolver,data);
%         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   FAIPAMAT 2.0 - 23/09/1999
% 
%   COMPUTES THE FEASIBLE DESCENT DIRECTION " d = d0 + rho*d1 "
%   when a Limited Memory Quasi - Newton Method is employed
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%disp('fdir_lim'),keyboard
fi=data(11);
alfa=data(12);
  
ntbind=neq+nbind+lenvlb+lenvub;
  
% Derivatives of the box constraints  
[gvlb,gvub]=boxderiv(nvar,lenvlb,lenvub,lvlb,lvub);

% Computing d0 and d1:

if isolver == 2   

% By solving the nonsymmetric system in (lambda,mu)
   
%   H=   MxHxM(eye(nvar),DELTA,GAMA,RMAT,DMAT,GAMAtGAMA,npairs)
      
   mat = MxHxM([gradg gvlb gvub],DELTA,GAMA,RMAT,DMAT,GAMAtGAMA,npairs);
   mat = diag([ones(neq,1);lambda])*mat;
   mat(neq+1:ntbind,neq+1:ntbind) =  mat(neq+1:ntbind,neq+1:ntbind) ...
                                    - diag(gt(neq+1:ntbind));
   vaux = Hxv(gradf,DELTA,GAMA,RMAT,DMAT,GAMAtGAMA,npairs);               
   mataux3 = mat\[(-diag([ones(neq,1);lambda])*[gradg gvlb gvub]'*vaux + ...
                  [gt(1:neq); zeros(ntbind-neq,1)]) ...
                  [omegaE; diag(lambda)*omegaI]];                                         
   
   lambda0 = mataux3(neq+1:ntbind,1);
   mu0 = mataux3(1:neq,1);           
   lambda1 = mataux3(neq+1:ntbind,2);
   mu1 = mataux3(1:neq,2); 
   
   glag0 =  (gradf + [gradg(:,neq+1:neq+nbind) gvlb gvub]*lambda0 ...
             + gradg(:,1:neq)*mu0);     
   d0 = - Hxv(glag0,DELTA,GAMA,RMAT,DMAT,GAMAtGAMA,npairs);
   glag1 = ([gradg(:,neq+1:neq+nbind) gvlb gvub]*lambda1 ...
              + gradg(:,1:neq)*mu1);       
   d1 = - Hxv(glag1,DELTA,GAMA,RMAT,DMAT,GAMAtGAMA,npairs);
          
elseif isolver == 3                         
   
% By solving the symmetric system in (lambda,mu)

   mat = MxHxM([gradg gvlb gvub],DELTA,GAMA,RMAT,DMAT,GAMAtGAMA,npairs);   
   mat(neq+1:ntbind,neq+1:ntbind) =  mat(neq+1:ntbind,neq+1:ntbind) ...
                                    - (diag(gt(neq+1:ntbind)./lambda));
   
   vaux = Hxv(gradf,DELTA,GAMA,RMAT,DMAT,GAMAtGAMA,npairs);               
   mataux3 = mat\[(- [gradg gvlb gvub]'*vaux + ...
                  [gt(1:neq); zeros(ntbind-neq,1)]) ...
                  [omegaE; omegaI]];                  
              
   lambda0 = mataux3(neq+1:ntbind,1);
   mu0 = mataux3(1:neq,1);           
   lambda1 = mataux3(neq+1:ntbind,2);
   mu1 = mataux3(1:neq,2); 
   
   glag0 =  (gradf + [gradg(:,neq+1:neq+nbind) gvlb gvub]*lambda0 ...
             + gradg(:,1:neq)*mu0);     
   d0 = - Hxv(glag0,DELTA,GAMA,RMAT,DMAT,GAMAtGAMA,npairs);
   glag1 = ([gradg(:,neq+1:neq+nbind) gvlb gvub]*lambda1 ...
              + gradg(:,1:neq)*mu1);       
   d1 = - Hxv(glag1,DELTA,GAMA,RMAT,DMAT,GAMAtGAMA,npairs);             
end   

% Updating c and the penalty function and its gradient

%   c = max(c,-1.2*mu0); MODIFICADO EM 4/09/2003 
   c = max(c,abs(1.2*mu0));
pen = f;
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
        
% Lagrangian's gradient norm

if isolver==3
   mglag0=norm(glag0);                     
else
   mglag0=norm(gradf+[gradg gvlb gvub]*[mu0;lambda0]); 
end

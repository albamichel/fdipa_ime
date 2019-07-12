%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function[d,d0,d02,lambda0,mu0,c,pen,gpen,mglag0] = fdir_box(B,f,gt,...
         gradf,gradg,lvlb,lvub,lambda,omegaI,omegaE,c,neq,nvar,...
         nbind,lenvlb,lenvub,isolver,errcg1,errcg2,data);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   FAIPAMAT 2.0 - 23/09/1999
% 
%
% COMPUTES THE FEASIBLE DESCENT DIRECTION " d = d0 + rho*d1 "

fi=data(11);
alfa=data(12);
  
ntbind=neq+nbind+lenvlb+lenvub;
nbox= lenvlb+lenvub;

% Calculus of the diag. matrix and the vector 
% corresponding to the contribution
% of the box constraints: 
%

[matbox,vecbox]=redbox(nvar,nbind,lambda,omegaI,...
                       gt(neq+nbind+1:neq+nbind+lenvlb),...
                  gt(neq+nbind+lenvlb+1:ntbind),lenvlb,lenvub,lvlb,lvub);
              
B = B - diag(matbox); % ATTENTION ! B is modified here but B is not in
                      %             the list of outputs of the function.

% Computing d0 and d1:
if isolver == 1
   % By solving the system in (d,lambda,mu):
    mataux1 = [B gradg ; ...  
              gradg(:,1:neq)' zeros(neq,ntbind-nbox) ; ...
              diag(lambda(1:nbind))*gradg(:,neq+1:neq+nbind)'...
              zeros(nbind, neq) diag(gt(neq+1:ntbind-nbox))] ...
             \[[-gradf; -gt(1:neq) ; zeros(ntbind-neq-nbox,1)] ...
              [ vecbox;-omegaE;-diag(lambda(1:nbind))*omegaI(1:nbind)]];
              
   lambda0 = mataux1(nvar+neq+1:nvar+neq+nbind,1);
   mu0 = mataux1(nvar+1:nvar+neq,1); 
   mu1 = mataux1(nvar+1:nvar+neq,2); 
   d0 = mataux1(1:nvar,1);
   d1 = mataux1(1:nvar,2);
   [lambda0]=lam_box(nvar,nbind,d0,lambda0,lambda,...
                     gt(neq+nbind+1:neq+nbind+lenvlb),...
                     gt(neq+nbind+lenvlb+1:ntbind),...
                     lvlb,lvub);
                  


   clear mataux1

elseif isolver == 2   
   % By solving the nonsymmetric system in (lambda,mu)
         
   mataux2 = B\gradg;
        
   mataux3 = [(diag(lambda(1:nbind))*gradg(:,neq+1:neq+nbind)'*...
               mataux2(:,neq+1:neq+nbind) -diag(gt(neq+1:neq+nbind))) ...
               diag(lambda(1:nbind))*gradg(:,neq+1:neq+nbind)'*...
               mataux2(:,1:neq); ...
               gradg(:,1:neq)'*mataux2(:,neq+1:neq+nbind) ...
               gradg(:,1:neq)'*mataux2(:,1:neq)]...
              \[[-diag(lambda(1:nbind))*mataux2(:,neq+1:neq+nbind)'*gradf; ...
                (-mataux2(:,1:neq)'*gradf + gt(1:neq))] ...
             [diag(lambda(1:nbind))*mataux2(:,neq+1:neq+nbind)'*vecbox + ...
                diag(lambda(1:nbind))*omegaI(1:nbind); mataux2(:,1:neq)'*vecbox+omegaE]];
          
   lambda0 = mataux3(1:nbind,1);
   mu0 = mataux3(nbind+1:neq+nbind,1);          
   lambda1 = mataux3(1:nbind,2);
   mu1 = mataux3(nbind+1:neq+nbind,2); 
   clear mataux2    
   clear mataux3
   d0 = - B\(gradf + gradg(:,neq+1:neq+nbind)*lambda0 ...
             + gradg(:,1:neq)*mu0);
   d1 = - B\(-vecbox + gradg(:,neq+1:neq+nbind)*lambda1 ...
             + gradg(:,1:neq)*mu1);
   [lambda0]=lam_box(nvar,nbind,d0,lambda0,lambda,...
                     gt(neq+nbind+1:neq+nbind+lenvlb),...
                     gt(neq+nbind+lenvlb+1:ntbind),...
                     lvlb,lvub);
                  

     
elseif isolver == 3                         
   % By solving the symmetric system in (lambda,mu)
   
       
   mataux2 = B\gradg;
           
   mataux3 = [(gradg(:,neq+1:neq+nbind)'*mataux2(:,neq+1:neq+nbind) ...
              -diag(gt(neq+1:neq+nbind)./lambda(1:nbind))) ...
               gradg(:,neq+1:neq+nbind)'*mataux2(:,1:neq); ...
               gradg(:,1:neq)'*mataux2(:,neq+1:neq+nbind) ...
               gradg(:,1:neq)'*mataux2(:,1:neq)]...
              \[[-mataux2(:,neq+1:neq+nbind)'*gradf; ...
               (-mataux2(:,1:neq)'*gradf + gt(1:neq))] ...
               [mataux2(:,neq+1:neq+nbind)'*vecbox + ...
               omegaI(1:nbind); mataux2(:,1:neq)'*vecbox+omegaE]];
          
   lambda0 = mataux3(1:nbind,1);
   mu0 = mataux3(nbind+1:neq+nbind,1);          
   lambda1 = mataux3(1:nbind,2);
   mu1 = mataux3(nbind+1:neq+nbind,2); 
   clear mataux2    
   clear mataux3
   d0 = - B\(gradf + gradg(:,neq+1:neq+nbind)*lambda0 ...
             + gradg(:,1:neq)*mu0);
   d1 = - B\(-vecbox + gradg(:,neq+1:neq+nbind)*lambda1 ...
             + gradg(:,1:neq)*mu1);
   [lambda0]=lam_box(nvar,nbind,d0,lambda0,lambda,...
                     gt(neq+nbind+1:neq+nbind+lenvlb),...
                     gt(neq+nbind+lenvlb+1:ntbind),...
                     lvlb,lvub);  
   
elseif isolver == 4                         
   % By eliminating lambda and solving a system in d0 and mu
    
   mataux3 = [(B - gradg(:,neq+1:neq+nbind)*...                                               
               diag(lambda(1:nbind)./gt(neq+1:neq+nbind))*...
               gradg(:,neq+1:neq+nbind)') ...
               gradg(:,1:neq); gradg(:,1:neq)' zeros(neq,neq)]...
              \[[-gradf; -gt(1:neq)] ...
               [(gradg(:,neq+1:neq+nbind)*...                                               
               diag(lambda(1:nbind)./gt(neq+1:neq+nbind))*omegaI(1:nbind)...
               + vecbox); -omegaE]];
                         
   d0 = mataux3(1:nvar,1);
   d1 = mataux3(1:nvar,2);
   mu0 = mataux3(nvar+1:nvar+neq,1);                                                    
   mu1 = mataux3(nvar+1:nvar+neq,2);
   lambda0 = - diag(lambda(1:nbind)./gt(neq+1:neq+nbind))*...
      gradg(:,neq+1:neq+nbind)'*d0;
   [lambda0]=lam_box(nvar,nbind,d0,lambda0,lambda,...
                     gt(neq+nbind+1:neq+nbind+lenvlb),...
                     gt(neq+nbind+lenvlb+1:ntbind),lvlb,lvub);
                                      
   clear mataux3   
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
        
% Lagrangian's gradient norm

glag0=gradf+gradg*[mu0;lambda0(1:nbind)]; 
j=0;
jj=0;
for i=1:nvar
   if lvlb(i)==1 
      j=j+1;
      glag0(i)= glag0(i)-lambda0(nbind+j);
   end
   if lvub(i)==1
      jj=jj+1;
      glag0(i)= glag0(i)+lambda0(nbind+lenvlb+jj);
   end
end                     
mglag0=norm(glag0);  % Modificado el 14/12/99 (de glag para glag0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function  [B,iterB]=updath(iterB,B,x,oldx,t,auxvec,gradf, ...
                           gradg,glag,rho2,glag2,bind,bindold, ...
                           mu0,lambda0,modg,nvar,ncstr,neq);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   FAIPAMAT 2.0 - 23/09/1999
% 
%   Updating of the Quasi - Newton matrix
%
%   B is an approximation of the inverse of the Hessian of the Lagrangian 
%
%   by exacting computing the inverse of BFGS  updating rule
% 
%   (see Luenberger pag. 269 eq. 40), modified by Powell to give
%
%   positive definite matrices,(Luenberger pag.60).
%
%   B computed here is exctly the inverse of B computed by updatb.

iterB = iterB+1;
if iterB == 3*nvar+1  
   B = eye(nvar,nvar);
   iterB=0;
end

delta = x-oldx;
gama=gradf+auxvec;
if neq >0
   for i=1:neq
      gama=gama + (gradg(:,i)/modg(i))*mu0(i);
   end
end
j=neq;
k=neq;
if ncstr > neq
   for i=1:ncstr-neq
      if bindold(i)==1, j=j+1;end
      if bind(i)==1, k=k+1;end
      if bindold(i)==1 & bind(i)==1
         gama=gama+(gradg(:,k)/modg(j))*lambda0(j-neq);
      end
   end
end




% Employs the inverse of B to apply Powell's modification

delgam = delta'*gama;
invBdel = -(t*glag + (t^2)* rho2 * glag2);
delinvBdel = delta'*invBdel;

if delgam < 0.2*delinvBdel 
   fifi = 0.8*delinvBdel/(delinvBdel - delgam);
   gama = fifi*gama + (1-fifi)*invBdel; 
   delgam = delta'*gama;
end
 
Bgam=B*gama;
gamBgam = gama'*Bgam; 
  
% Computes the dual of the BFGS matrix

%BB=B-(Bgam*Bgam')/(gamBgam)+(delta*delta')/delgam

% Computes the inverse of BFGS matrix

B = B + (1 + gamBgam/delgam)*(delta*delta'/delgam) - ...
        (delta*Bgam'+Bgam*delta')/delgam; 


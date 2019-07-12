%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function  [B,iterB]=updatb(iterB,B,x,oldx,auxvec,gradf, ...
                           gradg,bind,bindold,lambda0,mu0,...
                           modg,nvar,ncstr,neq);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   FAIPAMAT 2.0 - 23/09/1999
% 
%   Updating of the Quasi - Newton matrix
%
%   B is an approximation of the Hessian of the lagrangian by BFGS
%
%   updating rule (Luenberger pag. 269 eq. 36), modified by Powell to give
%
%   positive definite matrices,(Luenberger pag.60).


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


delgam = delta'*gama;  
Bdel=B*delta;
delBdel = delta'*Bdel;

if delgam < 0.2*delBdel 
   fifi = 0.8*delBdel/(delBdel - delgam);
   gama = fifi*gama + (1-fifi)*Bdel; 
   delgam = delta'*gama;
end

B=B-(Bdel*Bdel')/(delBdel)+(gama*gama')/delgam;

%delgam,eig(B),pause
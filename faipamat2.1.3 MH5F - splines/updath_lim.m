%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                 
function  [DELTA,GAMA,DMAT,RMAT,GAMAtGAMA,npairs]=...
          updath_lim(DELTA,GAMA,DMAT,RMAT,GAMAtGAMA,x,oldx,t,auxvec,gradf, ...
          gradg,glag,rho2,glag2,bind,bindold,mu0,lambda0,modg,nvar, ...
          ncstr,neq,npairs,mnpairs);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   FAIPAMAT 2.0 - 23/09/1999
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Updates the inverse of the BFGS Matrix by Limited Memory method
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
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
  
% Computes the inverse of the BFGS matrix by limited memory technique      

[DELTA,GAMA, DMAT,GAMAtGAMA,RMAT,npairs] = ...
            BFGS_lim(DELTA,GAMA, DMAT,GAMAtGAMA,RMAT,delta,gama,npairs,mnpairs);

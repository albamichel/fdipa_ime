%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function[f,g,bind,nbind,gradf,gradg,counter]=...
...
compute(x,g,fun,gfun,neq,neq_op,ncstr,icall,ilsearch,ifilt,bind,nbind,...
        counter, nprob,data,idata,iutil,rutil)
%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   FAIPAMAT 2.0 - 24/09/1999
% 
%
%
%  Controls computation of functions and derivatives
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% keyboard
nvar=length(x);
ideriv   = idata(14);
ifpoint  = idata(15);
indf_aux=0;
indgradf_aux=0;
if icall == 1, 

% 1) Computes f, all g (g_all)
% 2) CALLS FILTER 
% 3) Computes gradf, derivatives of equalities and binding 
%    inequality constraints.
%
% INPUT: x
%
% OUTPUT: f(x), g_all(x), gradf(x), gradg(1:neq+nbind)   

   indf=1;
   indg=ones(ncstr,1);
     
   indgradf=1;
   indgradg=ones(neq,1); 
   
   if ifilt==1 
      [bind,nbind]=filtcons(g,neq,ncstr);
      indgradg = [indgradg; bind]; 
   elseif ifilt==0 
      nbind = ncstr-neq;
      bind=ones(ncstr-neq,1);
      indgradg = [indgradg; ones(ncstr-neq,1)];
   end
    
   counter(1)=counter(1) + 1;
   counter(2)=counter(2) + 1;
   counter(3)=counter(3) + ncstr;
   counter(4)=counter(4) + neq + nbind;

elseif icall == 2 | icall == 4 | icall == 6, 

% Computes f, all the constraints: g_all, gradf and, if ilsearch=1,
% the derivatives of allthe equality constraints, to test the line 
% search criteria.
%
% INPUT: x
%
% OUTPUT: f(x), g_all(x), gradf(x), [gradg(x)(1:neq)] 

	gradf=zeros(nvar,0);gradg=zeros(nvar,0);

   indf=1;
   indg=ones(ncstr,1);
   indgradf=0;
   if ilsearch==1, indgradf=1; end
   indgradg = zeros(ncstr,1);
   if ilsearch==1, indgradg = [ones(neq,1); zeros(ncstr-neq,1)];end
   
   counter(1)=counter(1) + 1;
   if ilsearch==1,counter(2)=counter(2) + 1;end
   counter(3)=counter(3) + ncstr;
   if ilsearch==1,counter(4)=counter(4) + neq;end
   

elseif icall == 3,
   
% Computes the equality constraints and the binding inequality, to
% calculate d2
%
% 
%
% INPUT: x 
%
% OUTPUT: g1(x)(1:neq+nbind) 

	gradf=zeros(nvar,0);gradg=zeros(nvar,0);

   indf=0;
   indgradf=0;
   indgradg=zeros(ncstr,1);
       
   if ifilt==1 
      indg =[ones(neq,1);bind];
   elseif ifilt==0 
      nbind= ncstr-neq;
      indg = ones(ncstr,1);
   end
   
   counter(3)=counter(3) + neq + nbind;
     
elseif icall==5 % 

% 1) Filters 
% 2) Computes derivatives of the binding 
%    inequality constraints, to calculate a new d
%
%    The derivatives of the equality constraints remain 
%    unchanged
%
% 3) If ilsearch = 2 or 3, also computes the derivatives of the
%    objective and of the equality constraints.
%
% INPUT: x, 
%
% OUTPUT: if ilsearch=1: gradg(neq+1:neq+nbind)
%         if ilsearch=2 | 3: gradg(1:neq+nbind)


   f=[];gradf=zeros(nvar,0);gradg=zeros(nvar,0);
   indf=0;
   indg=zeros(ncstr,1);
   indgradf=0;

   if ifilt==1 
      [bind,nbind]=filtcons(g,neq,ncstr);
      indgradg = [zeros(neq,1); bind]; 
   elseif ifilt==0 
      nbind = ncstr -neq;
      bind=ones(ncstr-neq,1);
      indgradg = [zeros(neq,1); bind];
   end
         
   if ilsearch == 2 | ilsearch == 3,
      indgradf=1;
      indgradg(1:neq) = ones(neq,1);
   end   
  
   counter(4)=counter(4) + nbind;   
   if ilsearch == 2 | ilsearch == 3,
      counter(2)=counter(2) + 1;
      counter(4)=counter(4) + neq; 
   end
elseif icall == 10, 

% 1) Computes inequality constraints for the feasibility test
%
% INPUT: x
%
% OUTPUT: g_all(neq+1:ncstr)  

   indf=0;
   indg=[zeros(neq,1);ones(ncstr-neq,1)];
      
   f=0;gradf=zeros(nvar,0);gradg=zeros(nvar,0);   
   counter(3)=counter(3) + ncstr-neq;
 
elseif icall == 11, 

% 1) Computes f, all g(1:neq) (g_all)
% 2) CALLS FILTER 
% 3) Computes gradf, derivatives of equalities and binding 
%    inequality constraints.
%
% INPUT: x, g_all(neq+1:ncstr)
%
% OUTPUT: f(x), g_all(x), gradf(x), gradg(1:neq+nbind)   

   indf=1;
   indg=[ones(neq,1);zeros(ncstr-neq,1)];
   indgradf=1;
   indgradg=ones(neq,1); 
   
   if ifilt==1 
      [bind,nbind]=filtcons(g,neq,ncstr);
      indgradg = [indgradg; bind]; 
   elseif ifilt==0 
      nbind = ncstr-neq;
      bind=ones(ncstr-neq,1);
      indgradg = [indgradg; ones(ncstr-neq,1)];
   end
  
   counter(1)=counter(1) + 1;
   counter(2)=counter(2) + 1;
   counter(3)=counter(3) + neq;
   counter(4)=counter(4) + neq + nbind;
end
% If we are solving the aux. problem: 
if ifpoint~=0 & icall~=10
   indg = [zeros(neq_op,1);indg];
   indgradg = [zeros(neq_op,1);indgradg];
   indf_aux=indf; indgradf_aux=indgradf;
   if indf==1,indf=0; counter(1)=counter(1)-1;end
   if indgradf==1,indgradf=0; counter(2)=counter(2)-1;end
   z=x(nvar);
   x=x(1:nvar-1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if icall ~= 5
   [f,g]=feval(fun,x,indf,indg,nprob,iutil,rutil);
end 

if icall~=3 & icall~=10
   if ideriv==0
      [gradf,gradg]=feval(gfun,fun,x,indgradf,indgradg,nprob,iutil,rutil);
   else
      [gradf,gradg]=feval('findif',fun,x,indgradf,indgradg,nprob, ...
                          data,idata,iutil,rutil);
   end
end
% If we are solving the aux. problem:
if ifpoint~=0 & icall~=10 
   if indf_aux == 1
      f = z;
   end
   if indgradf_aux == 1
      gradf = [zeros(nvar-1,1); 1];
   end
   g=g-z;
   gradg(nvar,:) = -1;
   x=[x;z];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

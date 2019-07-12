%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [df,dg] = findif(fun,x,indgradf,indgradg,nprob,data,idata,iutil,rutil)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   FAIPAMAT 2.0 - 24/09/1999
%
%   Evaluates derivatives using finite differences
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  size(dg)= (nvar,ncstrder)
%
%  indgradf=0 - grad is not computed ( df=[] )
%  indgradf=1 - grad is computed
%
%  indgradg(i)=0 - dg(:,i) is not computed
%  indgradg(i)=1 - dg(:,i) is computed
%
%  length(indgradg)= ncstr
%
%  if indgradg(i)=0 for i=1:ncstr, take dg=zeros(nvar,0)

ideriv = idata(14);
h = data(8);

% ideriv ==1 => Central differences
% ideriv ==2 => Forward differences

% ncstrder = number of constraint derivatives

nvar=length(x); %number of variables
ncstr = length(indgradg);
ncstrder = indgradg'*ones(ncstr,1);
df=[];

if ncstrder > 0
   dg = zeros(0,ncstrder);
else
   dg=zeros(nvar,0);
end

if ideriv ==1, % Central Differences
   for i=1:nvar
      x(i)=x(i)+h;
      [ff,gf]=feval(fun,x,indgradf,indgradg,nprob,iutil,rutil);
		x(i)=x(i)-h-h;
      [fb,gb]=feval(fun,x,indgradf,indgradg,nprob,iutil,rutil);
      x(i)=x(i)+h;
      if indgradf
         df(i,1)=(ff-fb)/(h+h);
      end

      if ncstrder > 0
         dg_i=(gf-gb)'/(h+h);
         dg=[dg;dg_i];
      end
   end
elseif ideriv == 2, % Forward Differences

   [fx,gx]=feval(fun,x,indgradf,indgradg,nprob,iutil,rutil);
   for i=1:nvar
      x(i)=x(i)+h;
      [ff,gf]=feval(fun,x,indgradf,indgradg,nprob,iutil,rutil);
		x(i)=x(i)-h;
      if indgradf
         df(i,1)=(ff-fx)/h;
      end

      if ncstrder > 0
         dg_i=(gf-gx)'/h;
         dg=[dg;dg_i];
      end
   end
end


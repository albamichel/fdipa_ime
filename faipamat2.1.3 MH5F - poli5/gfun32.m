function[df,dg]=gfun32(fun,x,indgradf,indgradg,nprob,iutil,rutil)


%
%  size(dg)= (nvar,ncstr)
%
%  indgradf=0 - df is not computed ( df=[])
%  indgradf=1 - df is computed
%
%  indgradg(i)=0 - dg(:,i) is not computed
%  indgradg(i)=1 - dg(:,i) is computed
%
%  length(indgradg)= ncstr
%
%  if indgradg(i)=0 for i=1:ncstr, take dg=zeros(nvar,0)

df=[]; dg=zeros(3,0);
if indgradf==1
	df(1)=10*x(1)-2*x(2)+2*x(3);      
	df(2)=-2*x(1)+26*x(2)+6*x(3);     
	df(3)=2*(x(1)+3*x(2)+x(3));
end

j=0;
if indgradg(1)==1
	j=j+1;   
	dg(1,j)=1.0;
	dg(2,j)=1.0;
   dg(3,j)=1.0;
end   
if indgradg(2)==1
	j=j+1;   
	dg(1,j)=3*x(1)^2;
	dg(2,j)=-6;
   dg(3,j)=-4;
end

df=df(:);


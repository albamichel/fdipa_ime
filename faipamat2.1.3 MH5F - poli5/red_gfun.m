%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function[df,dg]=red_gfun(df,dg,indgradf,indgradg)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   FAIPAMAT 2.0 - 27/09/1999
%

nvar=size(df,1);
ncstr=size(dg,2);

if indgradf==0
   df=[];
end   

j=0;
for i=1:ncstr
   if indgradg(i)==1;
      j=j+1;
      dg(:,j)=dg(:,i);
   end   
end
if j < ncstr
	dg(:,j+1:ncstr)=[]; %zeros(nvar,0);   
end

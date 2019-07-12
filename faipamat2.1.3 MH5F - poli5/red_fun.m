%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function[f,g]=red_fun(f,g,indf,indg)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   FAIPAMAT 2.0 - 27/09/1999
%


ncstr=length(g);

if indf==0
   f=[];
end   

j=0;
for i=1:ncstr
   if indg(i)==1;
      j=j+1;
      g(j)=g(i);
   end   
end
if j < ncstr
   g(j+1:ncstr)=[];
end   
g=g(:);

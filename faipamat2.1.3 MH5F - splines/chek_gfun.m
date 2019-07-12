function[]=chek_gfun(gfun,x,indgradf,indgradg,gradf,gradg);


if indgradf==1 & isempty(gradf)==1
   disp(' ')
   disp(['ERROR IN ' gfun ])
   error('When "indgradf=1", the objective gradient "gradf" must be computed')
end
if indgradf==0 & isempty(gradf)==0
   disp(' ')
   disp(['ERROR IN ' gfun ])
   error('When "indgradf=0", if must be taken the objective gradient "gradf=zeros(length(x),0)"')
end
if size(gradg,2) ~= sum(indgradg)
   disp(' ')
   disp(['ERROR IN ' gfun ])
   error('The number of columns of "gradg" is iconsistent with "indgradg"')
end   
if indgradf==1 & length(gradf)~=length(x)
   disp(' ')
   disp(['ERROR IN ' gfun ])
   error('The length of "gradf" is not equal to "length(x)"')
end
if size(gradg,1)~= length(x)
   disp(' ')
   disp(['ERROR IN ' gfun ])
   error('The number of lines of "gradg" is not equal to "length(x)"')
end   

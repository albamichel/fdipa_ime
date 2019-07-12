function[]=chek_fun(fun,indf,indg,f,g);

if indf==1 & isempty(f)==1
   disp(' ')
   disp(['ERROR IN ' fun ])
   error('When "indf=1", the objective "f" must be computed')
end
if indf==0 & isempty(f)==0
   disp(' ')
   disp(['ERROR IN ' fun ])
   error('When "indf=0", if must be taken the objective "f=[]"')
end
if length(g) ~= sum(indg)
   disp(' ')
   disp(['ERROR IN ' fun ])
   error('The length o the contraints vector "g" is iconsistent with "indg"')
end   


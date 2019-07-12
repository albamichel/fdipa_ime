
function[f,g]=fun32(x,indf,indg,nprob,iutil,rutil)

%  fun32 - computes f and g 
% 
%  f is the objective function
%  g(ncstr) is a vector with the response contraints ordered as follows:
% 
%  g(1),...,g(neqlin)     - Linear equality constraints
%  g(neqlin+1),...,g(neq) - Nonlinear equality constraints
%  g(neq+1),...,g(ncstr)  - Inequality constraints (it does not
%                           include box constraints)
%
%  indf=0 - f is not computed ( Take f=[])
%  indf=1 - f is computed
%
%  indg(i)=0 - g(i) is not computed
%  indg(i)=1 - g(i) is computed
%
%  length(indg)= ncstr
%
%  if indg(i)=0 for i=1:ncstr, take g=[]
%
%
%  
%  iutil - Integer utility vector, employed to store data in
%                  probXX, funXX and gfunXX.
% 
%  rutil - Real utility vector, employed to 
%                  store data in probXX, funXX and gfunXX.
%
     
f=[];g=[];

if indf==1
   f=(x(1)+3*x(2)+x(3))^2+4*(x(1)-x(2))^2;
end

j=0;
if indg(1)==1
	j=j+1;   
	g(j)=-1+x(1)+x(2)+x(3);
end

if indg(2)==1
	j=j+1; 
	g(j)=x(1)^3-6*x(2)-4*x(3)+3;
end

g=g(:);

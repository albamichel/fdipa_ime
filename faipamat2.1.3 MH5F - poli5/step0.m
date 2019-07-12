%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [t,tbox]=step0(oldx,vlb,vub,lvlb,lvub,lambda0,lambda,tar, ...
                        g,oldgdg,d,d02,nvar,neq,lenvlb,lenvub,nbind,data,idata)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   FAIPAMAT 2.1 - 27/12/1999
%
%   Initial step for the line search in the Feasible Direction Method
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


eta  = data(15);

%   Box step

tbox = 1;
if lenvlb > 0
   for i =1:nvar
      if lvlb(i)==1 & d(i)<0 
         if (vlb(i)-oldx(i))==0
            d(i)=0;
         else 
            tvlb=(1-tar*eta)*(vlb(i)-oldx(i))/d(i);
            tbox = min(tbox, tvlb); 
         end
      end
   end 
end

if lenvub > 0
   for i=1:nvar
      if lvub(i)==1 & d(i) > 0 
         if(vub(i)-oldx(i))==0
            d(i) = 0; 
         else
            tvub=(1-tar*eta)*(vub(i)-oldx(i))/d(i);
            tbox = min(tbox, tvub);
         end
      end
   end 
end

t = tbox; 

if idata(5)==1,
   for  i = 1:nbind
      if lambda0(i)>0
         tlamb=(1-tar*eta)*lambda(i)/lambda0(i);
         t = min(t, tlamb);
      end 
   end 
end
   
if idata(6)==1,
   for i = neq+1:neq+nbind
      if g(i)<0           
         if oldgdg(i)>0
            tg=-(1-tar*eta)*g(i)/oldgdg(i);
            t = min(t, tg);
         end
      end        
   end 
end 

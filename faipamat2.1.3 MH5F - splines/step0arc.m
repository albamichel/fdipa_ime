%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [t,tbox]=step0arc2_0(oldx,vlb,vub,lvlb,lvub,lambda0,lambda,tar, ...
                        oldg,d,d2,d02,nvar,neq,lenvlb,lenvub,nbind,data,idata);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%   FAIPAMAT 2.0 - 27/09/1999
%
%   Initial step for the line search in the Feasible Arc Method

%   Box step

eta  = data(15);

tbox = 1;
if lenvlb > 0
   for i =1:nvar
      if lvlb(i)==1
         if (vlb(i)-oldx(i))==0
            d(i)=max(d(i),0);
            d2(i)=0;
         else
            tvlb=(1-tar*eta)* ...
                  minpos(roots([d2(i); d(i); oldx(i)-vlb(i)]));
            if length(tvlb) == 1    
               tbox = min(tbox, tvlb);
            end      
         end
      end
   end 
end
if lenvub > 0
   for i=1:nvar
      if lvub(i)==1  
         if(vub(i)-oldx(i))==0
            d(i)=min(d(i),0)
            d2(i) = 0; 
         else
            tvub=(1-tar*eta)*...
                  minpos(roots([d2(i); d(i); oldx(i)-vub(i)]));
            if length(tvub) == 1     
               tbox = min(tbox, tvub);
            end                          
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

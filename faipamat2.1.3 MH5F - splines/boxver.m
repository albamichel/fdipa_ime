%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
function[x,ibox] = boxver(x,nvar,vub,vlb,lvlb,lvub)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   FAIPAMAT 2.0 - 24/09/1999
%
%   VERIFIES INPUT DATA CORRESPONDING TO BOX CONSTRAINTS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    ibox=0;
    for i=1:nvar
     if lvlb(i) == 1 & lvub(i) == 1,
        if (vlb(i) > vub(i)), 
           disp('ERROR: WRONG STATEMENT OF THE BOUNDS');ibox=-1;return
        end
     end

     if lvlb(i)==1, 
       if ~(x(i)>vlb(i)),        
         disp('INFEASIBLE LOWER BOUND ');
         if lvub(i)==1,             
            x(i)=vlb(i)+abs(vub(i)-vlb(i))/100.0;
            disp(['WARNING: It was taken x0(', num2str(i),')=',num2str(x(i))])
         else
            delta=0.2*abs(vlb(i));
            if delta==0, delta=0.001; end
            x(i)=vlb(i)+delta;
            disp(['WARNING: It was taken x0(', num2str(i),')=',num2str(x(i))])
         end
       end 
     end
     if lvub(i)==1, 
         if ~(x(i)<vub(i)),
           disp('INFEASIBLE UPPER BOUND');
           if lvlb(i)==1,
               x(i)=vub(i)-abs(vub(i)-vlb(i))/100.0;
               disp(['WARNING: It was taken x0(', num2str(i),')=',num2str(x(i))])
           else
               delta=0.2*abs(vub(i));
               if delta==0, delta=0.001; end
               x(i)=vub(i)-delta;
               disp(['WARNING: It was taken x0(', num2str(i),')=',num2str(x(i))])
           end 
         end 
      end
    end
    
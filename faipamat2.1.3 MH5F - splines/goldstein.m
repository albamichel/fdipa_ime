%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [t,tl,tr,ftl,gtl_all,ftr,gtr_all,istalin,cod]=...
         goldstein(t,tbox,tl,tr,f,bind,g_all,oldf,...
                   oldg_all,oldgdf,oldgdg,ftl,gtl_all,...
                   ftr,gtr_all,tar,neqlin,neq,ncstr,data);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
%
%
%   FAIPAMAT 2.1 - 23/12/1999
%
%   GOLDSTEIN'S LINE SEARCH 
%
%   
%   WITH FILTERING          
%      
%

 

istalin=0;
told=t;                   
tbind=[ones(neq,1);bind]; % index of equality + binding inequality constr.

% PARAMETERS FOR STOPPING CRITERIA IN LINE SERACH: 
              

eta  = data(15);
eta1 = data(16);
eta2 = data(17); 
deps = 1.e-7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


if f <= oldf+t*eta1*oldgdf & (ncstr==neqlin | all(g_all(neqlin+1:ncstr)<=0))  
   
   cod = 1; % Verifies the upper test criterium
      
   if f >= oldf+t*eta2*oldgdf | ncstr==neqlin | ...
            any(g_all(neqlin+1:ncstr) >= eta*oldg_all(neqlin+1:ncstr)) ...
          | any(g_all(neqlin+1:ncstr) > -deps) | t==tbox | t==1 

      cod = 2; % Verifies the lower test criterium
      istalin = 1; % GOLDSTEINS'S CRITERIUM IS VERIFIED
   else
      cod = 3;      % EXTRAPOLATION                
      if tr > 0
         cod = 4;   % EXTRAPOLATION WITH tr > 0
         tnew = tr;
         tf = tl + inquaf3P(ftl,f,ftr,t-tl,tr-tl); 
        
         if tf <= 2*t, tf = 2*t + .2*(tr-t); end % Do not accept tf < 2*t
        
         tnew = min(tnew,tf);        
         if ncstr > neqlin
            for i=neqlin+1:ncstr 
               tnew = min(tnew,tl + inquag3P(gtl_all(i),g_all(i),gtr_all(i),...
                          tar*eta*oldg_all(i),t-tl,tr-tl));
            end
         end
      else  
         cod = 5;  % EXTRAPOLATION WITH tr = 0
         tnew = inf;                      
         tf=inquadf2P(f,oldf,oldgdf,t);
         if tf <=1.1*t, tf = 1.1*t; end
         tnew = min(tnew,tf);
         j=neqlin;                      
         if ncstr > neqlin
            for i=neqlin+1:ncstr
               if tbind(i)==1
                  j=j+1;
                  tg=inquag2P(g_all(i),oldg_all(i),...
                              tar*eta*oldg_all(i),oldgdg(j),t);
               else
                  tg=inlin(g_all(i),oldg_all(i),...
                           tar*eta*oldg_all(i),t);
               end
               if tg<=t, tg=tnew; end
               tnew = min(tnew,tg);
            end                                  
         end
      end     
      t = min([tnew, tbox,1]);
      if abs(t-told)< .2*told,      
         t=1.2*told;
         t=min([t,tbox,1]);
      end 
      if (tr>0 & abs(t-tr)< 1.e-4),
         t=0.5*(tl+tr);
      end            
      tl = told; ftl = f; gtl_all = g_all; 
   end   
else % (f <= oldf+t*eta1*oldgdf & all(g_all(neqlin+1:ncstr)<=0))

   cod = 6; % INTERPOLATION
   tnew = t - tl; % Modified 08/03/99, 

   if f > (oldf + t*eta1*oldgdf)
      if tl ==0,
          tf=inquadf2P(f,oldf,oldgdf,t);
      else
          tf = tl + inquaf3P(ftl,f,ftr,t-tl,tr-tl); 
      end   
      tnew = min(tnew,tf); 
   end
   if ncstr > neqlin
      j=0;  
      for i=1:ncstr
         if tbind(i)==1
            j=j+1;
         end  
         if g_all(i) >= 0
            if tl == 0
               if tbind(i)==1
                  tg=inquag2P(g_all(i),oldg_all(i),...
                              tar*eta*oldg_all(i),oldgdg(j),t);                    
               else
                  tg=inlin(g_all(i),oldg_all(i),...
                           tar*eta*oldg_all(i),t);
               end
            else
               tg=inquag3P(gtl_all(i),g_all(i),gtr_all(i),...
                   tar*eta*oldg_all(i),t-tl,tr-tl);
            end     
            tnew = min(tnew,tg);
         end %(if g_all(i) > 0)
      end   %(for i=1:ncstr)
   end     % (if ncstr > 0)
   tr = told; ftr = f; gtr_all = g_all;
   tnew=tnew+tl;
  
   t = min([tnew, tbox,1]);
   if abs(told-t)< 1.e-4*(told-tl),      
      t=tl+.9*(told-tl);
   end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [t,tl,tr,ftl,gtl_all,gdftl,ftr,gtr_all,gdftr,istalin,cod]=...
         wolfe(t,tbox,tl,tr,f,bind,g_all,gdf,oldf,...
                oldg_all,oldgdf,oldgdg,ftl,gtl_all,gdftl,...
                ftr,gtr_all,gdftr,tar,neqlin,neq,ncstr,data);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%
%
%   FAIPAMAT 2.1 - 27/12/1999
%
%   WOLFE'S LINE SEARCH <<< WITHOUT >>> CONSTRAINTS DERIVATIVES
%   
%   WITH FILTERING          
%      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
istalin=0;
told=t;                   
tbind=[ones(neq,1);bind]; % index of equality + binding inequality constr.

% PARAMETERS FOR STOPPING CRITERIA IN LINE SERACH: 
              
eta  = data(15);
eta1 = data(16);
eta2 = data(17); 
deps = 1.e-7   ;

if f <= oldf+t*eta1*oldgdf & (ncstr==neqlin | ...
                             all(g_all(neqlin+1:ncstr)<=0)) 
   
   cod = 1; % Verifies the upper test criterium
      
   if gdf >= eta2*oldgdf | ncstr==neqlin | any(g_all(neqlin+1:ncstr) ...
                           >= eta*oldg_all(neqlin+1:ncstr)) ...
          | any(g_all(neqlin+1:ncstr) > -deps) | t==tbox | t==1 

      cod = 2;     % Verifies the lower test criterium
      istalin = 1; % WOLFE'S CRITERIUM IS VERIFIED
   else
      cod = 3;      % EXTRAPOLATION                
      if tr > 0
         cod = 4;   % EXTRAPOLATION WITH tr > 0          
         tnew = tr;
         tnew = min(tnew,t+incubf(ftr,f,gdftr,gdf,tr-t)); 
         if ncstr > neqlin %>>>>>>
            for i=neqlin+1:ncstr 
                         
               tnew = min(tnew,tl+inquag3P(gtl_all(i),g_all(i),gtr_all(i),...
                          tar*eta*oldg_all(i),t-tl,tr-tl));
                          
            end
         end
      else  
         cod = 5;  % EXTRAPOLATION WITH tr = 0
         tnew = inf;                      
         tf=incubf(f,oldf,gdf,oldgdf,t);
         tnew = min(tnew,tf);
         j=neqlin;                      
         if ncstr > neqlin %>>>>
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
%      if abs(t-told)< 1.e-1,   % 9/10/2000
%keyboard
%         t=2.*told; 
%        if (tr>0 & abs(t-tr)< 1.e-4),
%            t=0.5*(tl+tr);
%         end
%         t=min([t,tbox,1]);
%      end             
      tl = told; ftl = f; gtl_all = g_all; gdftl = gdf;  
   end   
else % (f <= oldf+t*eta1*oldgdf & all(g_all(neqlin+1:ncstr)<=0))
   cod = 6; % INTERPOLATION
   
   %   tnew = inf;
   tnew = t - tl; 
   if f > (oldf + t*eta1*oldgdf)     
      tf=incubf(f,ftl,gdf,gdftl,t-tl);
      tnew = min(tnew,tf); 
   end
   if ncstr > neqlin %>>>>>> 19/7/2000
      j=0;  
      for i=neqlin+1:ncstr  %>>>>>>>>17/4/2001
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

                            %if t-tl > .05*(tr-tl) % Modified 22/12/1999
               tg=inquag3P(gtl_all(i),g_all(i),gtr_all(i),...
                           tar*eta*oldg_all(i),t-tl,tr-tl);
                            %else
                            %tg=inlin(g_all(i),gtr_all(i),...
                                      % tar*eta*oldg_all(i),tr);    
                           % end       
            end     
            tnew = min(tnew,tg);
         end %(if g_all(i) > 0)
      end   %(for i=1:ncstr)
   end     % (if ncstr > 0)
   tr = told; ftr = f; gtr_all = g_all;gdftr = gdf; 
   tnew=tnew+tl;
   t = min([tnew, tbox,1]);
   if abs(told-t)< 1.e-1*(told-tl),   
      t=tl+.5*(told-tl);
   end 
end

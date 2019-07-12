
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

function [x,f,g,lambda0,mu0,counter] = ...
...
faipa(fun,gfun,x,vlb,vub,nvar,ncstr,neq,neqlin,lvlb,lvub,nprob,...
      data,idata,iutil,rutil);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
%
% FAIPA - FEASIBLE ARC INTERIOR POINT ALGORITHM 
%
%   FAIPAMAT 2.1 - 23/12/1999
% 
%   Makes "FILTERING' of the inequality constraints
%
%  "Normalizes" the constraints, at each iteration
%
%   Reduces the box constraints
%
%   Solves the internal linear systems:
%
%   1 - The primal-dual system
%   2 - The dual system
%   3 - The dual symmetric system
%   4 - The primal system + dual variables corresponding
%       to the equality constraints
%   5 - The dual symmetrized with the conjugate gradient
%       method
%
%   Quasi - Newton Updating Rule: BFGS (Luenberger, pag ...)
%
%   Approximates the Hessian of the Lagrangian when:
%
%   1 - The primal-dual system is solved
%   2 - The primal + dual (inequality) system is solved
%   3 - In all the cases when the box constraints are reduced
%
%   Approximates the inverse when:
%
%   1 - The dual system is solved
%   2 - The dual symmetric system is solved
%
%   Limited memory Quasi - Newton method
%   Aproximates the inverse with the same rule than before, when:
%
%   1 - The primal-dual system is solved
%   2 - The primal  + dual (inequality) system is solved
%   Only when the box constraints are not reduced
%   
%
%   Line search:
%
%   1 - Wolfe
%   2 - Goldstein
%   3 - Armijo
%
%  All vectors are column
%
%  gradg - Derivative matrix of the constraints. 
%          The derivative of each constraint is in a column.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         

disp(' ######  FAIPAMAT2.1.1 #####')

imethod  = idata(1);
ifilt    = idata(2);
iscale   = idata(3);
iupdateb = idata(4);
isolver  = idata(10);
iredbox  = idata(11);
ilsearch = idata(12);
ideriv   = idata(14);
ifpoint  = idata(15);
 
errcg1   = data(6);
errcg2   = data(7);
delfd    = data(8);

if imethod==1, disp('FEASIBLE DIRECTION INTERIOR POINT ALGORITHM'), end
if imethod==2, disp('FEASIBLE ARC INTERIOR POINT ALGORITHM'), end

if iscale==0, disp('THE CONSTRAINTS ARE NOT SCALED'),end
if iscale==1, disp('THE CONSTRAINTS ARE SCALED'),end

if ifilt==0, disp('THE CONSTRAINTS ARE NOT FILTERED'),end
if ifilt==1, disp('THE CONSTRAINTS ARE FILTERED'),end

if iupdateb==0, disp('First Order Algorithm'),end
if iupdateb==1, disp('Quasi-Newton Algorithm'),end
if iupdateb==2, 
   disp('Limited Memory Quasi-Newton Algorithm'),
   mnpairs = idata(13);
   disp(sprintf('Maximum number of stored pairs delta-gama = %3.0f', mnpairs))
   if isolver ~= 2 & isolver ~= 3
      error('Please take isolver =2 or isolver =3')
   end
end

if isolver==1 | isolver==4 | iredbox==1           
   disp('B approximates the Hessian of the Lagrangian')
else            
   disp('B approximates the inverse of the Hessian of the Lagrangian')
end   

disp('INTERNAL LINEAR SYSTEMS:')
if     isolver==1, disp('Solves a system in (d,lambda,mu)')
elseif isolver==2, disp('Solves a nonsymmetric system in (lambda,mu)')
elseif isolver==3, disp('Solves a symmetric system in (lambda,mu)')
elseif isolver==4, disp('Solves a system in d0 and mu')
elseif isolver==5, 
  disp('Stopping Criteria for the CG:')
  disp(sprintf(' errcg1 = %0.1g; errcg2 = %1.1g', errcg1,errcg2))
else   disp('Please, define isolver = 1, 2, 3,4 or 5'), return
end
if isolver == 5 & iredbox == 1
   error('Do not take isolver = 5 and iredbox =1')
end   

if iredbox==1, disp('THE BOX CONSTRAINTS ARE REDUCED'),end
if iredbox==0, disp('THE BOX CONSTRAINTS ARE NOT REDUCED'),end
            % If iupdateb = 2,TAKE:  iredbox = 0 !!
if iredbox==1 & iupdateb==2
   error('If iupdateb = 2, TAKE: iredbox = 0')
end
if ilsearch==1, disp('Wolfe Criteirum in the line search'),end
if ilsearch==2, disp('Goldstein Criterium in the line search'),end
if ilsearch==3, disp('Armijo Line Search  '),end 

if ideriv == 1, disp('Derivation by Central Finite Diferences'), end
if ideriv == 2, disp('Derivation by Forward Finite Diferences'), end
if ideriv==1 | ideriv==2, 
   disp(sprintf('Increment for Finite Diferences: delfd =  %0.5g', delfd))
end

% Actualizar en FORTRAN90  >>> 12/07/2000
if ncstr==neq & ifpoint==1
   disp('The problem has no inequality constraints. Takes ifpoint=0.')
   idata(15)=0;ifpoint=idata(15);
end 
if ncstr==neq & ifpoint==2
   error('The problem has no inequality constraints. Take ifpoint=0')
end 
% Actualizar en FORTRAN90  >>> 12/07/2000   
if ifpoint ==1, 
   disp('Makes a search of an initial feasible point and optimizes')
end

if ifpoint ==2, 
   disp('Makes a search of an initial feasible point and stops')
end

% nbind - NUMBER OF INEQULITY BINDING CONSTRAINTS

initial=1;
% initial
%            1 - initial point
%            0 - not initial point

istatus=0;
% istatus
%            0 - initialization
%            1 - feasible direction d + initial step 
%            2 - feasible arc + initial step
%            3 - newstep

istop=0;
counter=[0 0 0 0 0]';
% if istop~=0, stops iterations 

% main iterates counter:
iter=0;
counter(5)=0;


% Counter of functions and derivatives calculus
%
% counter(1)=number of evaluations of the objective function
% counter(2)=number of evaluations of the objective function's derivatives
% counter(3)=number of constraints evaluations (each const. counted)
% counter(4)=number of evaluations of the constraints' derivatives
%            (each const. counted)
% counter(5)=Number of iterations, including the search for an interior
%            point.

while (istop==0 | istop == 10) %>>>>>>> 16/12

   if istatus==0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                %
%    INITIALIZATION                                              %                                                            %
%                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      neq_op = neq; %Stores "neq" of the original problem %>>16/12/1999
      neqlin_op = neqlin;
      if ncstr==neq & ifilt==1 
         ifilt=0;
         disp('Since there are not inequality constraints, takes ifilt=0')
      end
      
      rho2=0; d2=zeros(nvar,1);glag2=zeros(nvar,1);
      

      if iscale==1 & ncstr==0,
         iscale=0;
         disp('ncstr=0 - Takes iscale=0')
      end 
      
      glow=[];
      gup=[];

      istalin=0;

      [x,ibox] = boxver(x,nvar,vub,vlb,lvlb,lvub);
      
      if ibox==-1,return,end
      
% COUNTS THE BOX CONSTRAINTS

      lenvlb=0;
      lenvub=0;

      for i=1:nvar
         if lvlb(i)==1
            lenvlb = lenvlb+1;
         end
         if lvub(i)==1
            lenvub = lenvub+1;
         end
      end
         
      ntcstr=ncstr+lenvlb+lenvub;
      nbox  = lenvlb+lenvub;      
                       

% INITIAL B

      B=[];   
      if iupdateb == 1 | iupdateb == 0
         B = eye(nvar,nvar);
      end

% Counter of B's updatings, after reinitialization

      if iupdateb==1, 
         iterB=0;
      end 

% Initialitation for Limited Memory

      npairs=0; DDELTA=zeros(nvar,0); GGAMA=zeros(nvar,0);... %16/12/1999
      DMAT=[]; RMAT=[]; GAMAtGAMA=[]; 



% INITIAL OMEGA

      omegaE=zeros(0,1);
      omegaI_all=zeros(0,1);
      omegaI=zeros(0,1);
      
      if neq>0, omegaE = ones(neq,1); end
      if neqlin>0, omegaE(1:neqlin,1) = zeros(neqlin,1); end
      if ntcstr > neq, omegaI_all = ones(ntcstr-neq,1); end
      

% INITIAL PENALTY PARAMETER c
 
      c = zeros(neq,1);
      
      istatus=1; % Goes to the calculus of the feasible direction
      
      g=[];bind=[];nbind=[];
      
%     CALL 1



      if ifpoint==0 %%>>> 16/12

         icall=1; %%>>> 16/12


         [f,g_all,bind,nbind,gradf,gradg,counter]=...
         ...
         compute(x,g,fun,gfun,neq,neq_op,ncstr,icall,ilsearch,ifilt,bind,nbind,...
                 counter, nprob,data,idata,iutil,rutil);
      else
      
         icall = 10;

         if istop==0
      
            [f,g_all,bind,nbind,gradf,gradg,counter]=...
            ...
            compute(x,g,fun,gfun,neq,neq_op,ncstr,icall,ilsearch,ifilt,...
                    bind,nbind,counter, nprob,data,idata,iutil,rutil);
         end 

         if any(g_all(1:ncstr-neq)>=-eps)
            disp('The initial point is infeasible')
            disp('Searching a Feasible Point') % 20/7/2000 >>>>   			
         
            [x,g_all,lvlb,lvub,nvar,ncstr,ntcstr,neq,neqlin,B,DDELTA,GGAMA,...
                         d2,glag2,omegaE,c] = ...
            ...
            set_aux_prob(x,g_all,lvlb,lvub,nvar,ncstr,ntcstr,neq,neqlin,...
                         B,DDELTA,GGAMA,d2,glag2,omegaE,c,idata);
                         %>>>>>>> 30/08/2000
       
         else
            
            idata(15)=0;ifpoint=idata(15);
            disp(' >>>>> The given initial point is feasible')
            
         end
         
         icall = 11;
         
         g_all=[zeros(neq,1);g_all];
         [f,g_all(1:neq),bind,nbind,gradf,gradg,counter]=...
         ...
         compute(x,g_all(1:neq),fun,gfun,neq,neq_op,ncstr,icall,ilsearch,ifilt,bind,nbind,...
                 counter,nprob,data,idata,iutil,rutil);                      
                                              

      end      
      
   elseif istatus==1     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
%                                                                 %%%%
%     Calculus of the Feasible Direction                          %%%%
%                                                                 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


         iter=iter+1;
         counter(5)=counter(5)+1;
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
%                                                                    %
%   FOR THE INTIAL POINT  ONLY                                       %
%                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                     %							
         if initial==1                                               %                                                                                                           %
                                                                     %                                                                     %         
% FEASIBILITY TEST                                                   % 
                                                                     %	  
            if any(g_all(neq+1:ncstr)>=0),                           %   
               disp('INITIAL POINT INFEASIBLE'); return              %
            end                                                      %
                                                                     %
            head(nprob,nvar,ncstr,neq,nbox,f,fun)                    % 
                                                                     %
                                                                     %	  
% SIGN OF THE EQUALITY CONSTRAINTS                                   %
%                                                                    %	  
% To make negative the equality constraints at the initial point     %
                                                                     %	  
                                                                     %	  
           if neq>0                                                  %
              chsig = -sign(g_all(1:neq));                           %
              for i=1:neq                                            %
                 if chsig(i)==0, chsig(i)=1; end                     %
              end                                                    %  
           end                                                       %
                                                                     %	
% Changes the sign of the equality constraints                       %
                                                                     %
            if neq>0                                                 %
               g_all(1:neq)=diag(chsig)*g_all(1:neq);                %
               gradg(:,1:neq)=gradg(:,1:neq)*diag(chsig);            % 
            end                                                      %
         end                                                         %	
                                                                     %	         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


         
% UPDATING QUASI-NEWTON MATRIX B 
                 
         if iupdateb==1 & initial==0
            if isolver==1 | isolver==4 | iredbox==1
            
               % B approximates the Hessian
            
               [B,iterB]=updatb(iterB,B,x,oldx,gradf,gradfold, ...
                                gradg,gradgold,bind,bindold,mu0,lambda0,...
                                modg,nvar,ncstr,neq);
            else
            
                % B approximates the inverse of the Hessian
                                                
                [B,iterB]=updath(iterB,B,x,oldx,t,gradf,gradfold, ...
                                gradg,gradgold,glag,rho2,glag2,bind,bindold,mu0,lambda0,...
                                modg,nvar,ncstr,neq);
                                               
            end                                               
         end
         if iupdateb==2 & initial==0
            
            % Limited memory approximation of the inverse of the Hessian

            [DDELTA,GGAMA,DMAT,RMAT,GAMAtGAMA,npairs]=...
               updath_lim(DDELTA,GGAMA,DMAT,RMAT,GAMAtGAMA,x,oldx,t,gradf,gradfold, ...
               gradg,gradgold,glag,rho2,glag2,bind,bindold,mu0,...
               lambda0,modg,nvar,ncstr,neq,npairs,mnpairs);           
         end   

% MOUNTS THE BINDING CONSTRAINTS VECTOR g(x)

         j=0;g=g_all(1:neq);
         for i=1:ncstr-neq
            if bind(i)==1;
               j=j+1;
               indexbin(j)=neq+i;  % Index of binding inequality constraints
               g=[g;g_all(neq+i)];
            end    
         end
         g=g(:);  
                         
      
% SCALINF OF THE BINDING CONSTRAINTS
         
% Computes the module of the binding derivatives

         modg=[];
         if iscale==0
            modg=ones(neq+nbind,1);
         end   
         if iscale==1 & neq > 0
            for i=1:neq
               modg(i)=norm(gradg(:,i));
               if modg(i)<=1.e-10, modg(i)=1.e-10; end
            end 
         end
         if iscale==1 & ncstr > neq                      
            j=neq;
            for i=1:ncstr-neq
               if bind(i)==1
                  j=j+1;
                  modg(j)=norm(gradg(:,j));
                  if modg(j)<=1.e-10, modg(j)=1.e-10; end 
               end    
            end
         end   
            
% Scales binding constraints

         if iscale==1 & neq > 0,
            for i=1:neq
               g(i)=g(i)/modg(i);
               g_all(i)=g(i);
               gradg(:,i)=gradg(:,i)./modg(i);                
            end
         end
         if iscale==1 & nbind > 0,                      
            for i=neq+1:neq+nbind
               g(i)=g(i)/modg(i);
               g_all(indexbin(i-neq))=g(i);
               gradg(:,i)=gradg(:,i)./modg(i); 
            end    
         end 
         
% Mounts binding omegaI
         omegaI=zeros(0,1);
         j=0;
         for i=neq+1:ncstr
            if bind(i-neq)==1
               omegaI=[omegaI; omegaI_all(i-neq)];
            end    
         end

         omegaI = [omegaI; omegaI_all(ncstr-neq+1:ntcstr-neq)];
      
      
% BOX CONSTRAINTS CALCULUS

         glow=zeros(0,1); gup=zeros(0,1);
         j=0;
         k=0;
         for i=1:nvar
            if lvlb(i)==1
               j=j+1;    
               glow(j) = vlb(i)-x(i);
            end
            if lvub(i)==1
              k=k+1; 
              gup(k) = x(i)-vub(i);
            end
         end 
               
% TOTAL BINDING CONSTRAINTS VECTOR 
      
         gt=[g;glow(:);gup(:)];
         
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
%                                                                         %
%   FOR THE INTIAL POINT  ONLY                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                          %	                        
                                                                          %      
% INITIAL LAMBDA_ALL                                                      %  
                                                                          %
         if initial==1 % For the initial point only                       %
            initial=0;                                                    %
            lambdaS=1e+5;                                                 %
            lambda_all = -1./min([g_all(neq+1:ncstr);glow'; gup'],...     %
                                 -1/lambdaS);                             %
            lambda_all=lambda_all(:);                                     %                          % 
         end                                                              %                                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
% MOUNTS THE BINDING LAMBDAS VECTOR 

         j=0;lambda=zeros(0,1);
         if ncstr>neq     
            for i=1:ncstr-neq
               if bind(i)==1;
                  lambda=[lambda; lambda_all(i)];
               end    
            end
         end   
         
         lambda = [lambda; lambda_all(ncstr-neq+1:ntcstr-neq)];
     
%         
%      
% COMPUTATION OF THE FEASIBLE DIRECTION " d" 

         if iredbox == 0
            if iupdateb == 0 | iupdateb == 1                    
                  [d,d0,d02,lambda0,mu0,c,pen,gpen,glag,mglag0] = ...
               fdir(B,f,gt,gradf,gradg,lvlb,lvub,lambda,omegaI,omegaE,c,...
                    neq,nvar,nbind,lenvlb,lenvub,isolver,errcg1,errcg2,data);                    
            else            
                  [d,d0,d02,lambda0,mu0,c,pen,gpen,glag,mglag0] = ...
               fdir_lim(DDELTA,GGAMA,RMAT,DMAT,GAMAtGAMA,npairs,f,gt,...
                        gradf,gradg,lvlb,lvub,lambda,omegaI,omegaE,c,neq,nvar,...
                        nbind,lenvlb,lenvub,isolver,data);                     
            end 
         else
         
               [d,d0,d02,lambda0,mu0,c,pen,gpen,mglag0] = ...
            fdir_box(B,f,gt,gradf,gradg,lvlb,lvub,lambda,omegaI,omegaE,c,...
                     neq,nvar,nbind,lenvlb,lenvub,isolver,errcg1,errcg2,data);
         end

         modd0=d02^.5;
         tar=min(.5, d02);
keyboard                                       
% SAVING OLD VALUES  (At the last main iterate, that is , for t=0)
  
         oldx = x; oldpen = pen; oldg_all = g_all; oldg=g; 
         oldgdg = gradg'*d; oldgdpen = gpen'*d; 

% INITIAL STEP

         if imethod ==1
         
            [t,tbox]=step0(oldx,vlb,vub,lvlb,lvub,lambda0,lambda,tar,...
                           g,oldgdg,d,d02,nvar,neq,lenvlb,lenvub,nbind,data,idata);                           
                                              
         else 
          
            t=1;
             
         end
       
     
% NEWPOINT

          x = oldx + t * d;         
                                     
          if imethod==1
          
             iterlin=0; % Iterations in the line search
             
             scod=[];
          
             istatus = 3; % Goes to the line search iterations
             
          else
                           
             istatus = 2; % Goes to the calculus of the centering direction
             
          end 
          
               
% Stores the old derivative and old binding

          gradgold = gradg;  % SCALED
          gradfold = gradf;
          bindold = bind;
          nbindold=nbind;                
          
          if imethod ==1 
           
             
% COMPUTING the objective, all the constraints
% If ilsearch=1 also computes the derivatives of the
% objective and of the equality constraints

%            CALL 2


	          icall=2;
	          g_all=[];

            [f,g_all,bind,nbind,gradf,gradg,counter]=...
            ...
            compute(x,g_all,fun,gfun,neq,neq_op,ncstr,icall,ilsearch,ifilt,bind,nbind,...
                    counter, nprob,data,idata,iutil,rutil);
          else
          
% CALL 3
%
% Computing the equality and binding inequality constraints to get "d2".

            icall=3;
      
            g1=[];
            
            [fdum,g1,bind,nbind,gradfdum,gradgdum,counter]=...
            ...
            compute(x,g1,fun,gfun,neq,neq_op,ncstr,icall,ilsearch,ifilt,bind,nbind,...
                    counter, nprob,data,idata,iutil,rutil);
         
          end          
                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
%                                                                 %%%%
%     Calculus of the Feasible Arc                                %%%%
%                                                                 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
      
   elseif istatus == 2
      
% Changes the sign of equality constraints
                  
         if neq>0            
            g1(1:neq)=diag(chsig)*g1(1:neq);            
         end 
      
% Scaling the constraints
      
         if iscale==1
            for i=1:neq
               g1(i)=g1(i)/modg(i);
            end   
            for i=1:nbind
               g1(neq+i)=g1(neq+i)/modg(neq+i);
            end     
         end
         
% Setting OMEGAE2, OMEGAI2

         omegaI2=zeros(0,1);         
 
         if nbind>0,
            omegaI2= g1(neq+1:neq+nbind) - oldg(neq+1:neq+nbind) ...
                  - oldgdg(neq+1:neq+nbind);
         end
          omegaE2=zeros(0,1);
  
         if neq>0          
            omegaE2= g1(1:neq) - oldg(1:neq) - oldgdg(1:neq);
            if neqlin>0
               omegaE2(1:neqlin,1) = zeros(neqlin,1); 
            end            
         end
    
% Makes equal to zero the negative elements of omega
%
%        for i=1:neq
%           if omegaE2(i)<0, omegaE2(i)=0; end
%        end
%
%        for i=1:nbind
%           if omegaI2(i)<0, omegaI2(i)=0; end  
%        end
       
% COMPUTATION OF THE CENTERING DIRECTION "d2"    

         if iredbox == 0
            if iupdateb == 0 | iupdateb == 1                   
               [d2,glag2] = cendir(B,gt,gradg,lvlb,lvub,lambda,omegaI2,...
                                   omegaE2,neq,nvar,nbind,lenvlb,lenvub,...
                                   isolver,errcg1,errcg2);
            else
               [d2,glag2] = cendir_lim(DDELTA,GGAMA,RMAT,DMAT,GAMAtGAMA,npairs,gt,...
                                       gradg,lvlb,lvub,lambda,omegaI2,omegaE2,...
                                       neq,nvar,nbind,lenvlb,lenvub,isolver);
            end                                               
         else                          
            [d2] = cendir_box(B,gt,gradg,lvlb,lvub,lambda,omegaI2,...
                              omegaE2,neq,nvar,nbind,lenvlb,lenvub,...
                              isolver,errcg1,errcg2);
         end                                               
                      
 % WE DON'T LET norm(d2) > norm(d0)
            
         modd2=norm(d2);
         rho2=0.;		
         if modd2 ~=0, rho2=min(1,modd0/modd2); end 
         d2=rho2*d2;
  
    
% INITIAL STEP ALONG  A FEASIBLE ARC


         [t,tbox]=step0arc(oldx,vlb,vub,lvlb,lvub,lambda0,lambda,tar, ...
                           g,d,d2,d02,nvar,neq,lenvlb,lenvub,nbind,data,idata);
                                                     
                

% COMPUTES A NEW POINT FOR THE INITIAL STEP                           
                           
                           
          x = oldx + t*d +t^2*d2;

          istatus=3; % Goes to the line search iterations
          
          iterlin=0; % Iterations in the line search counter
          
          scod=[];                                                 


% COMPUTING the objective, all the constraints
% If ilsearch=1 also computes the derivatives of the
% objective and of the equality constraints

% CALL 4

         icall=4;
         g_all=[];

         [f,g_all,bind,nbind,gradf,gradg,counter]=...
         ...
         compute(x,g_all,fun,gfun,neq,neq_op,ncstr,icall,ilsearch,ifilt,bind,nbind,...
                    counter, nprob,data,idata,iutil,rutil);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
%                                                                 %%%%
%     I)  VERIFIES LINE SEARCH CRITERIA                           %%%%
%                                                                 %%%%
%     II) COMPUTES A NEW STEP                                     %%%%
%                                                                 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
         
   elseif  istatus ==3
      

% Initial tr and tl for Wolfe's and Goldstein's line search        
      
         if iterlin==0         
            tl=0; pentl=oldpen; gtl_all=oldg_all; gdpentl=oldgdpen; gdgtl=oldgdg;
            tr=0; pentr=0; gtr_all=gtl_all; gdpentr=gdpentl; gdgtr=gdgtl;         
         end         
         
         if neq>0 
            g_all(1:neq) = diag(chsig)*g_all(1:neq);
            if ilsearch==1, gradg(:,1:neq) = gradg(:,1:neq)*diag(chsig);end
         end
         
% Maximum error in equality constraints

         erreq=0;
         if neq>0
            erreq=norm(g(1:neq),inf);
         end
         
% MOUNTS THE BINDING CONSTRAINTS VECTOR

         indexbin=[];
         g=g_all(1:neq);
         for i=1:ncstr-neq
            if bind(i)==1;
               indexbin=[indexbin; neq+i];  % Index of binding constraints
               g=[g;g_all(neq+i)];
            end    
         end
         g=g(:);            

% Scales the binding constraints at the new step

         if iscale==1 & neq > 0
            for i=1:neq
               g(i)=g(i)/modg(i);
               g_all(i)=g(i);      
            end
         end   

         if iscale==1 & nbind > 0,
            for i=neq+1:neq+nbind
               g(i)=g(i)/modg(i);
               g_all(indexbin(i-neq))=g(i);               
            end     
         end
     
% Scales the equality constraints derivatives at the new step, 
% (only for Wolfe's search)

         if iscale==1 & ilsearch==1 & neq >0,
            for i=1:neq          
               gradg(:,i)=gradg(:,i)./modg(i); 
            end     
         end         

% New penalty function in line search
    
         pen = f;
         if ilsearch==1, gpen = gradf; end
         if neq > 0
            pen  = pen - g(1:neq)'*c;
            if ilsearch==1, gpen = gpen - gradg(:,1:neq)*c; end
         end
         
% Derivatives with respect to the step-length t

         
         if imethod==1 & ilsearch == 1
            gdpen  = gpen'*d;
            gdg    = gradg'*d;
         end   
         if imethod==2 & ilsearch == 1  
            gdpen  = gpen'*(d + 2*t*d2);
            gdg    = gradg'*(d + 2*t*d2);
         end
         
% New step in line search
                    
         if  ilsearch == 1
                              
            [t,tl,tr,pentl,gtl_all,gdpentl,pentr,gtr_all,...
             gdpentr,istalin,cod]=...
            wolfe(t,tbox,tl,tr,pen,bind,g_all,gdpen,oldpen,...
                  oldg_all,oldgdpen,oldgdg,pentl,gtl_all,gdpentl,...
                  pentr,gtr_all,gdpentr,tar,neqlin,neq,ncstr,data);
           
         elseif  ilsearch == 2        

             [t,tl,tr,pentl,gtl_all,pentr,gtr_all,istalin,cod]=...
             goldstein(t,tbox,tl,tr,pen,bind,g_all,oldpen,...
                        oldg_all,oldgdpen,oldgdg,pentl,gtl_all,...
                        pentr,gtr_all,tar,neqlin,neq,ncstr,data);
         else                        

            [t,istalin,cod]=armijo(t,pen,g_all,oldpen,oldgdpen,neqlin,ncstr,data);
            
         end

  
                                                       
% String with the line search  history
                
         scod=[scod int2str(cod)]; 
                
         if istalin==1, % The line search criterium is verifyed.
         
% Un-scaling the binding constraints

            if iscale==1 & neq > 0,
               for i=1:neq
                  g(i)=g(i)*modg(i);
                  g_all(i)=g(i);
               end     
            end
            if iscale==1 & nbind > 0,
               for i=neq+1:neq+nbind
                  g(i)=g(i)*modg(i);
                  g_all(indexbin(i-neq))=g(i);
               end     
            end
            
% Un-scaling the equality constraints derivatives

            if iscale==1 & ilsearch==1 & neq > 0,
               for i=1:neq
                  gradg(:,i)=gradg(:,i).*modg(i); 
               end     
            end            
                          
            disp(sprintf('%3.0f %10.6g %10.6g %10.6g %10.6g %10.6g %4.0f %4s', ...
                          iter,f,modd0,mglag0,erreq,t,nbind,scod))
                       
  % TEST  OF STOPPING CRITERIA

            [istop]=stoptest(istop,modd0,mglag0,f,g_all,oldpen,pen,erreq,...
                             iter,iterlin,data,idata);
                             %%>>> 16/12
                                                              
            
            if istop==0,
                                    
               istatus=1; % Goes to the computation of a new feasible direction
            
% UPDATES            
            

     
% Updating the dual variables lambda
 
               [lambda]=updatlam(lambda0,d02);
               
% Mounting lambda_all % Analisar el posible scaling de lambda

               j=neq; lambda_all=[];
               for i=neq+1:ncstr
                  if bind(i-neq)==1
                  j=j+1;
                     lambda_all(i-neq)=lambda(j-neq);
                  else
                     lambda_all(i-neq) = -1./min(g_all(i),-1/lambdaS); % Dikin's updating
                  end
               end
               lambda_all=lambda_all(:);
               
               lambda_all = [lambda_all; lambda(nbind+1:nbind+nbox)];             
               
               
% Calling the filtering and the new derivatives of the binding
% inequality constraints.
% If ilsearch = 2 or 3, also computes the derivatives of the
% objective and of the equality constraints.

% CALL 5
               icall=5;
               
               bind=[]; nbind=0;
               if ilsearch == 1,
                  [fdum,gdum,bind,nbind,gradfdum,gradg(:,neq+1:neq+nbind),counter]=...
                  ...
                  compute(x,g_all,fun,gfun,neq,neq_op,ncstr,icall,ilsearch,ifilt,bind,nbind,...
                          counter, nprob,data,idata,iutil,rutil);
               end          
               
               if ilsearch == 2 | ilsearch ==3,
                  [fdum,gdum,bind,nbind,gradf,gradg,counter]=...
                  ...
                  compute(x,g_all,fun,gfun,neq,neq_op,ncstr,icall,ilsearch,ifilt,bind,nbind,...
                          counter, nprob,data,idata,iutil,rutil);
                   if neq>0 
                        gradg(:,1:neq) = gradg(:,1:neq)*diag(chsig);
                   end
         

               end       

                         
            else %(if istop==0) 
            
               dispstop(istop)  %%>>> 16/12
               if ifpoint==1
                  istatus=0;
                  initial=1;
                  iter=0;
                  [x,g_all,lvlb,lvub,nvar,ncstr,neq,neqlin] = ...
                  ...
                  reset_prob(x,g_all,lvlb,lvub,nvar,ncstr,neq,...
                             neqlin,neq_op,neqlin_op);
               end
            end  %(if istop==0)
                                      
         else   %(if istalin==1)                             
            if iterlin >= idata(9)
               error('ERROR, the maximum number of iterations in line search')
            else          
               iterlin=iterlin + 1;
            end 
            if imethod == 1,
               x = oldx + t*d;
            else
               x = oldx + t*d + t^2*d2;
            end
            
% CALL 6
 

% COMPUTING the objective, all the constraints
% If ilsearch=1 also computes the derivatives of the
% objective and of the equality constraints

            icall=6;

            g_all=[];
            
       	   [f,g_all,bind,nbind,gradf,gradg,counter]=...
       	   ...
        	   compute(x,g_all,fun,gfun,neq,neq_op,ncstr,icall,ilsearch,ifilt,bind,nbind,...
                     counter, nprob,data,idata,iutil,rutil);

         end   %(else if istalin==1)         
   end  % (else if istatus==3  )   
end % (while istop==0)  

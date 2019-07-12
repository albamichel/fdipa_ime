%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function head(nprob,nvar,ntcstr,neq,nbox,f,fun);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   disp(' ');
   disp(' ');
   disp(' ');
   disp('=============================================================')
   if isempty(fun)==1,
    disp(sprintf('      PROBLEMA TESTE No. %3.0f', nprob));
   else
    disp(fun);    
   end
   disp(sprintf('                        nvar=%3.0f  ncstr=%3.0f  neq=%3.0f  nbox=%3.0f', ...
         nvar,ntcstr,neq,nbox));
   disp('-------------------------------------------------------------')
   disp('ITER   OBJECTIVE   |d0|     |lagran|   erreq   STEP   nbind  icod');
   disp('=============================================================')
   disp(' '); 

   disp(sprintf('%3.0f %10.6g', 0, f));

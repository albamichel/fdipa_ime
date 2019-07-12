%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function  [auxvec]=update_aux(auxvec,gradg,bind,bindold,lambda0,ncstr,neq);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  Auxiliar vector for updatb (second part) 30/10/2000
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%keyboard
j=neq;
k=neq;
if ncstr > neq
   for i=1:ncstr-neq
      if bindold(i)==1, j=j+1;end
      if bindold(i)==1 & bind(i)==1
         auxvec=auxvec - gradg(:,j)*lambda0(j-neq);
      end
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function[matbox,vecbox]=redbox(nvar,nbind,lambda,omegaI,glow,gup,...
 lenvlb,lenvub,lvlb,lvub);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   FAIPAMAT 2.0 - 27/09/1999
%
%   When the box constraints are reduced:
%
%   Calculus of the diag. matrix corresponding to the contribution
%   of the box constraints: 
%
%   matbox=gradg*G^(-1)*diag(lambda)*gradg'
% 
%   Calculus of the vector corresponding to the contribution
%   of the box constraints: 
%
%   vecbox=gradg*G^(-1)*diag(lambda)*omegaI
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matbox=zeros(nvar,1);

j=0;
jj=0;
for i=1:nvar
   if lvlb(i)==1 
      j=j+1;
      matbox(i)= matbox(i)+1/glow(j)*lambda(nbind+j);
   end
   if lvub(i)==1
      jj=jj+1;
      matbox(i)= matbox(i)+1/gup(jj)*lambda(nbind+lenvlb+jj);
   end
end

vecbox=zeros(nvar,1);

j=0;
jj=0;
for i=1:nvar
   if lvlb(i)==1 
      j=j+1;
      vecbox(i)= vecbox(i)-1/glow(j)*lambda(nbind+j)*omegaI(nbind+j);
   end
   if lvub(i)==1
      jj=jj+1;
      vecbox(i)= vecbox(i)+...
         1/gup(jj)*lambda(nbind+lenvlb+jj)*omegaI(nbind+lenvlb+jj); % >>> ERROR corregido 
                                                                    % el 12/12/2000: gup()
   end
end
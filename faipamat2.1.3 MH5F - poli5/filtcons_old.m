%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function[bind,nbind]=filtcons(g,neq,ncstr)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   FAIPAMAT 2.0 - 24/09/1999
% 
%   DEFINES THE SET OF INEQUALITY BINDING CONSTRAINTS
%
%   This function is only a sample.
%
%   For any aplication the user can define the most
%   efficient criteria to select the binding constraints
%
%   bind - Is a vector of zeros and ones 
%   size(bind)-ncstr-neq
%   bind(i) = 1 => The (inequality) constraint g(neq+i) is 
%                  considered to compute the line serach

bind = zeros(ncstr-neq,1);
nbind=0;
for i=neq+1:ncstr
   if g(i) > -1.
      bind(i-neq)=1;
      nbind=nbind+1;
   end
end 
bind=bind(:);

%disp('filtcons')

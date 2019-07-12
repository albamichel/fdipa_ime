%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function[bind,nbind]=filtconsFP(z,g,ncstr,iutil,rutil)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   26/10/2000 
%   
% 
%   DEFINES THE SET OF INEQUALITY BINDING CONSTRAINTS FOR
%
%   THE AUXILIARY PROBLEM TO FIND A FEASIBLE POINT
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

bind = zeros(ncstr,1);
nbind=0;
for i=1:ncstr
   if g(i)+ z > -.1
      bind(i)=1;
      nbind=nbind+1;
   end
end 
bind=bind(:);

%disp('filtconsFP')

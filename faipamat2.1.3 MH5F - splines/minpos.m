%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function minpr = minposreal(v)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   FAIPAMAT 2.0 - 27/09/1999
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = length(v);

for i=1:n
   if imag(v(i))~=0, v(i)=0; end
end

v = max(v,0);
v = nonzeros(v);
minpr=min(v);  
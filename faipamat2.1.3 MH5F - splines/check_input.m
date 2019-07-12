%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function[] = check_input(x,vlb,vub,nvar,ncstr,neq,neqlin,lvlb,lvub)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Ckecks input data to FAIPA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
disp( 'Checking consistency of input data:')
disp('  ')
index=0;
if length(x) ~= nvar
   error( 'The length of x0  must be equal to nvar')
   index=1;
end   
if length(vlb) ~= nvar
   error( 'The length of vlb  must be equal to nvar')
   index=1;
end 
if length(vub) ~= nvar
   error( 'The length of vub  must be equal to nvar')
   index=1;
end 
if length(lvlb) ~= nvar
   error( 'The length of lvlb  must be equal to nvar')
   index=1;
end 
if length(lvub) ~= nvar
   error( 'The length of lvub  must be equal to nvar')
   index=1;
end 
if neq > ncstr
   error('It must be ncstr >= neq')
   index=1;
end
if index==0
   disp('"Input data is consistent"')
end   
disp('  ')
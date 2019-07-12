%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [S, Y, D, YtY, R1,l] = BFGS_lim(S, Y, D, YtY, R1,s, y,l,m)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   FAIPAMAT 2.1 -23/12/1999
% 
%******************************************************************************%
%     Updates the Hessian matrix implicitely defined according with the method %
%     L-BFGS                                                                   %
%                                                                              %
%     Nº of products: 2mn + mm + m                                             %
%                                                                              % 
%     s,y = Last pair of vectors                                               % 
%     l = Number of stored pairs of vectors s,y                                %
%     S, Y, D, YtY, R1 = Auxiliar matrices                                     %
%                => último par de vetores                                      %
%                                                                              %
%     Reference: Byrd, R. H., Nocedal, J. and Schnabel, R. B., "Representation %
%     of Quasi-Newton Matrices and Their Use in Limited Memory Methods",       %
%     Technical Report CU-CS-612-92, 1992, University of Colorado at Boulder   %
%******************************************************************************%

if l == m
    S   = S(:,2:m);
    Y   = Y(:,2:m);
    D   = D(:,2:m);
    YtY = YtY(2:m,2:m);
    R1=R1(2:m,2:m);
else
    l   = l + 1;
end

aux1 = y' * y;
aux2 = Y' * y;
YtY  = [YtY, aux2; aux2', aux1];
aux1 = s' * y; % aux1=1/rho
aux2 = -(1 / aux1) * (R1 * (S' * y));

R1   = [R1, aux2; zeros([1, (l - 1)]), (1 / aux1)];
D    = [D, aux1];
S    = [S, s];
Y    = [Y, y];



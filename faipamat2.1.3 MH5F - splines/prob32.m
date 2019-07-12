
%
%     Test Problem Nr. 32 - Schittkowski et al.
%

clear all;
format long;

nvar=3;    % Nr. of Variables
ncstr=2;   % Nr. of Constraints, excluded box contraints
neq=1;     % Nr. of Equality Constraints
neqlin=1;  % Nr. of Linear Equality Constraints
nprob=32;  % Nr. of the problem to be solved, required when
           % there are several problems in one file.

x0 = [.1 0.7 0.2 ]'; % Initial Point
vlb = [0 0 0];        % Lower Bound value in Box Constraints
vub = [999 999 999];  % Upper Bound value in Box Constraints
 
lvlb=[1 1 1];
lvub=[0 0 0];

% lvlb(i)=1 if there is a Lower Bound on x(i).
% lvlb(i)=0 if there is NOT a Lower Bound on x(i). Then,
%           vlb(i) is not taken in consideration.


% lvub(i)=1 if there is an Upper Bound on x(i).
% lvub(i)=0 if there is NOT an Upper Bound on x(i). Then,
%           vub(i) is not taken in consideration.

fun='fun32';   % File with the MATLAB function that computes
               % the Objective and the Constraints.

gfun='gfun32'; % File with the MATLAB function that computes
               % the Objective and the Constraints DERIVATIVES.

iutil=[]; % Integer utility vector, if needed by fun & gfun.
rutil=[]; % Real utility vector, if needed by fun & gfun.

[data,idata]=fdata; % CONTROL PARAMETERS

[x,f,g,lambda0,mu0,counter] = faipa(fun,gfun,x0,vlb,vub,nvar,ncstr,...
                                    neq,neqlin,lvlb,lvub,nprob,...
                                    data,idata,iutil,rutil)
% COMPUTED OPTIMUM:
%
% x ------- Variables
% f ------- Objective Function
% g ------- Constraints
% lambda0 - Lagrange Multipliers of Inequality Constraints.
% mu0 ----- Lagrange Multipliers of Equality Constraints

%
% counter(1)=number of evaluations of the objective function
% counter(2)=number of evaluations of the objective function's derivatives
% counter(3)=number of constraints evaluations (each const. counted)
% counter(4)=number of evaluations of the constraints' derivatives
% (each const. counted)
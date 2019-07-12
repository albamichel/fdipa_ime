    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [data,idata]=fdata( ); 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   FAIPAMAT 2.1 - 23/12/1999
% 
%   INPUT OF CONTROL DATA DEFINED BY THE USER
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Stopping criteria

data(1)=0.00001;       % Stopps if mlag < data(11)
data(2)=0.00001;       % Stopps if modd0   < data(12)
data(3)=0;             % Stopps if the decrease of f
                      % is less than data(13)
data(4)= 1.0e-5;      % Maximum accepted error for the equality 
                      % constraints
                   
data(5)= 0.00001;        % Stopps if "fun < data(5)"
                       % When looking for a feasible point take
                       % data(5).le.0

% Conjugate Gradient Solver

data(6)= 1.e-6;  %= errcg1; Stopping criterium for the Conjugate Gradient
data(7)= 1.e-9;  %= errcg2; Stopping criterium for the Conjugate Gradient

% Finite Differences Derivatives

data(8)= eps^(1/3); %= hfd; Increment for derivatives computation by 
                    % finite differences.
                    %
                    % Recommended values for h:
                    %                          1) Central FD: hfd = eps^(1/3)
                    %                          2) Forward FD: hfd = eps^(1/2)                
% Search Direction Calculus

data(11) = 1.;  % = phi   - Default: phi=1;
data(12) = .7;  % = alpha - Default: alpha=0.7;

% Line Search Parameters

data(15) = .5; % = eta  - Default = .5
data(16) = .1; % = eta1 - Default = .1
data(17) = .7; % = eta2 - Default = .7

data(18) = .8; % = nu   

% 
idata(1)=1;    % =imethod; = 1-FEASIBLE DIRECTION, =2 FEASIBLE ARC
idata(2)=0;    % =ifilt; =1-FILTERS, =0- DO NOT
idata(3)=0;    % =iscale; =1-SCALES CONSTRAINTS =0 DO NOT
idata(4)=1;    % =iupdateb; 
               % = 0 - Takes B = I
               % = 1 - Updates B by a quasi - Newton rule.
               % = 2 - Limited memory quasi - Newton update
               % If iupdateb = 2, PLEASE TAKE:
               %                         isolver = 2 or = 3 !!
               %                         iredbox = 0        !!

idata(5)=0;    % = 1 - CONTROL OF tlamb in the initial step
idata(6)=0;    % = 1 - CONTROL OF tg in the initial step


idata(7)=700;  % itmax = maximum number of iterations
idata(8)=4; %   ;  % itmin = minimum number of iterations
                   % Recommended: idata(8)=4
                   %              idata(8)=0 when looking for 
                   %              a feasible point
                   
idata(9)=400;  % itlinmax= maximum number of iterations 
               % in line search

idata(10)=1 ;  % =isolver; 
% Computes de Feasible Direction and/or the Feasible Arc:
%
% 1 - By solving the system in (d,lambda,mu)
% 2 - By solving the nonsymmetric system in (lambda,mu)
% 3 - By solving the symmetric system in (lambda,mu)
% 4 - By eliminating lambda and solving a system in d0 and mu
% 5 - By solving the symmetric system in (lambda,mu) with the
%     Conjugate Gradient Method
%
%     IF ISOLVER = 5, PLEASE TAKE IREDBOX =0 !!  

idata(11)=0; % = iredbox
             % = 1 - reduces box constraints
             % = 0 - do not

             % IF ISOLVER = 5, TAKE:  iredbox = 0 !!
             % If iupdateb = 2,TAKE:  iredbox = 0 !!
               
idata(12)=3; % ilsearch - defines line search:
             % = 1 - Wolfe's criteium 
             % = 2 - Goldsteins' criterium
             % = 3 - Armijo's line search            
             
idata(13)=6; % = mnpairs - Maximum number of pairs delta-gama 
             % stored for Limited Memory Quasi-Newton

idata(14)=1; % ideriv - Derivation
             % ideriv = 0, Derivates are computed by gfunXX.m
             % ideriv = 1, Central Finite Differences
             % ideriv = 2, Forward Finite Differences
             
idata(15)=1; % ifpoint =1: - Search for an initial feasible point
                             % and optimizes  
             % ifpoint =2: - Search for an initial feasible point
                             % and stops                                                          

% To find an initial interior point take "idata(5)=0" and,
% eventually, "data(4)=Inf" and data(5)=0
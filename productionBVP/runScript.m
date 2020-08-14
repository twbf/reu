%---------------RUNSCRIPT--------------------------%

addpath('BC_FVM/')
addpath('BCproj/')
addpath('HodoSolve/')

%%% Global variables:
global cf2 d xc FS IN LW
global a amp td b g g2 dx h N M
global H1 H2 c1 c2 x1 x2
global t0 Tf x0 Xf
global eta1 u1 Psi_nT Psi_nB
global phi psi numLam numSig

FS = 'FontSize';
IN = 'Interpreter';
LS = 'LineStyle';
LW = 'LineWidth';

%%% Physical parameters:
g  = 9.8;	% gravity acceleration
g2 = 0.5*g;	% g/2
cf2 = 0.0;	% friction coefficient

a  = -2.0;	% the left boundary (incident wave)
b  = 10.0;	% the right boundary (wall)
% Bottom slope:
td = 10.0/10.0;

%%% Initial condition parameters:
H1 = 0*0.006;
H2 = 0.018;
c1 = 0*0.4444;
c2 = 4.0;
x1 = 0*4.1209;
x2 = 1.6384;

% initial time and final time (domain of t)
t0 = 0;
Tf = 5;

% Hodograph Solver resolution parameters
numSig = 6;
numLam = 10000;

% domain of x
x0 = -2;
Xf = 10;

N  = 2000;					% number of grid points
M  = 1000;	% number of time instances where we project solution

[eta1, u1] = BCpull();   % solution via FVM
[Psi_nT, Psi_nB] = BCproj();   % hodohraph solutions via Data Projection
[phi, psi] = HodoSolve(); % hodograph Solver









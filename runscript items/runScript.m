%---------------RUNSCRIPT--------------------------%
% Runs:
% runwave.m
% HodographSolver.m
% Psi_n_Chebfun.m
%---------------------------------------------------%
clc;
close all;

global M N a b d td g g2 
global td H1 H2 c1 c2 x1 x2
global FS IN LS LW

M = 2000        % number of time steps
N = 2000        % number of spatial steps

g = 9.81;        % acceleration due to gravity
g2 = g*0.5;
td = 10.0/10.0;  % bathymetry slope

a = -2;           % spatial bounds in (x,t) plane
b = 10;
d = 1;             % char. depth

%%% Initial condition parameters:

H1 = 0.006;
H2 = 0;
c1 = 0.4444;
c2 = 0;
x1 = 4.16;
x2 = 0;


[u1_dim,eta1_dim] = run_wave(M,N,a,b,d,td,g,g2)





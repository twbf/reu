%notes: non-zero_velocity case

clear
close all
format longE

%%%libraries we use
addpath('FVM_numerical/');
addpath('NOAA_analytic/');
addpath('2018_analytic/');

%%% Global variables:
%all realted to initial conditions, bathymetry and time
global H1 H2 c1 c2 x1 x2 eta_0 eta_prime u_0 u_prime td t0 Tf x0 Xf upS lowS upL lowL K

%%% bathymetry and time

td = 10.0/10.0; %slope of bathymetry

t0 = 0; %initial time
Tf = 1.9; %final time (10 in Catalina 1 experiement)

x0 = -0.1;
Xf = 3;

%%%initial conditions

H1 = 0.006;
H2 = 0.018;
c1 = 0.4444;
c2 = 4.0;
x1 = 4.1209;
x2 = 1.6384;

%ranges are high just so there are no problems
eta_0 = chebfun(@(x) H1*exp(-c1*(x - x1).^2) - H2*exp(-c2*(x - x2).^2), [-2.5, 30]);
eta_prime = diff(eta_0);

%u_0   = chebfun(@(x) -2.*(sqrt(x+eta_0(x)+3)-sqrt(x+3)), [0 30]);
%u_prime = diff(u_0);

%resolution of solution
res = 1000;

%Nicolsky 2018 solution

% integration parameters
lowS = 0;  %lower bound on s
upS = Xf+1;   %upper bound on s
lowL = t0;  %lower bound on lambda
upL = 10;   %upper bound on lambda

K = 20; % upper bound on integration parameter

[ana_eta, ana_u] = CG_transform(true, 'psi_phi_projection_test4.mat', res, true, false);

[num_eta, num_u] = run_num();
%numeric solution

t0 = 0.01


disp(' ');
disp('Comparing');

%comparison
[diff, l2, ana_eta] = stat(num_eta, ana_eta, res, true);

disp(' ');
disp('Saving and Ploting');

%display
Plot_xt(num_eta, res, 1, '$\eta$ FVM', 'x', 't', 'eta');
Plot_xt(ana_eta, res, 2, '$\eta$ Nicolsky 2018', 'x', 't', 'u')
Plot_xt(diff, res, 3, '$\eta$ FVM - $\eta$ Nicolsky', 'x', 't', 'diff')


figure(4)
plot(linspace(t0, Tf, res), l2);
title(['$L2$ Norm Difference Between FVM and Nicolsky 2018'],'Interpreter','latex');
xlabel('$t$','Interpreter','latex');
ylabel('$L2$ Norm','Interpreter','latex');

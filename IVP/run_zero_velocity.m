
clear
close all
format longE

%%%libraries we use
addpath('FVM_numerical/');
addpath('NOAA_analytic/');
addpath('2018_analytic/');

%%% Global variables:
%all realted to initial conditions, bathymetry and time
global eta_0 eta_prime u_0 u_prime td t0 Tf x0 Xf upS lowS upL lowL K

%%% bathymetry and time

td = 10.0/10.0; %slope of bathymetry

t0 = 0; %initial time
Tf = 10; %final time (10 in Catalina 1 experiement)

x0 = 0;
Xf = 10;

%%% initial conditions

H1 = 0.006;
H2 = 0.018;
c1 = 0.4444;
c2 = 4.0;
x1 = 4.1209;
x2 = 1.6384;

eta_0 = chebfun(@(x) H1*exp(-c1*(x - x1).^2) - H2*exp(-c2*(x - x2).^2), [-3, 50]);
eta_prime = diff(eta_0);

u_0   = chebfun(@(x) 0, [-3 21]);
u_prime = diff(u_0);

%resolution of solution
res = 1000;


%%% catalina 1 solution
%[cat_eta, cat_u] = catalina_transform(true, 1000);

%%% numeric solution
%[num_eta, num_u] = run_num();

%catalina 1 solution
[cat_eta, cat_u] = catalina_transform('catalina1_phi_psi2.mat', res);

%%% Nicolsky 2018 solution
% integration parameters
lowS = 0.0001;  %lower bound on s
upS = Xf+1;   %upper bound on s
lowL = t0;  %lower bound on lambda
upL = 10;   %upper bound on lambda

K = 30; % upper bound on integration parameter

[ana_eta, ana_u] = CG_transform(true,'zero_initial_u_catalina_k30.mat', res, false, true);

t0 = 0.0001; %initial time
Tf = 1.9; %final time (10 in Catalina 1 experiement)

x0 = -0.1;
Xf = 3;


%%% numeric solution
[num_eta, num_u] = run_num();

%comparison
[d, l2, ana_eta] = stat(num_eta, ana_eta, res, true); %nuclear_option

[d_cat, l2_cat, cat_eta] = stat(num_eta, cat_eta, res, true); %nuclear_option

%display
Plot_xt(num_eta, res, 1, '$\eta$ FVM', 'x', 't', '$\eta$')
Plot_xt(cat_eta, res, 2, '$\eta$ NOAA Catalina1', 'x', 't', '$\eta$')
Plot_xt(ana_eta, res, 3, '$\eta$ Nicolsky', 'x', 't', 'u')

Plot_xt(d, res, 4, '$\eta$ FVM - $\eta$ Nicolsky', 'x', 't', 'diff')
Plot_xt(d_cat, res, 5, '$\eta$ FVM - $\eta$ NOAA Catalina1', 'x', 't', 'difference')

figure(6)
plot(linspace(t0, Tf, res), l2);
title(['$L2$ Norm Difference Between FVM and Nicolsky 2018'],'Interpreter','latex');
xlabel('$t$','Interpreter','latex');
ylabel('$L2$ Norm','Interpreter','latex');

figure(7)
plot(linspace(t0, Tf, res), l2_cat);
title(['$L2$ Norm Difference Between FVM and NOAA Catalina 1'],'Interpreter','latex');
xlabel('$t$','Interpreter','latex');
ylabel('$L2$ Norm','Interpreter','latex');

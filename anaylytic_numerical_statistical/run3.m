%notes: non-zero_velocity case

clear
close all
format longE

%%%libraries we use
addpath('numerical/');
addpath('catalina1/');
addpath('anaylytic/');

%%% Global variables:
%all realted to initial conditions, bathymetry and time
global H1 H2 c1 c2 x1 x2 eta_0 eta_prime u_0 u_prime td t0 Tf x0 Xf

%%% bathymetry and time

td = 10.0/10.0; %slope of bathymetry

t0 = 0; %initial time
Tf = 1; %final time (10 in Catalina 1 experiement)

x0 = 0;
Xf = 10;

%%%initial conditions

H1 = 0.006;
H2 = 0.018;
c1 = 0.4444;
c2 = 4.0;
x1 = 4.1209;
x2 = 1.6384;

eta_0 = @(x) H1*exp(-c1*(x - x1).^2) - H2*exp(-c2*(x - x2).^2);
eta_prime = chebfun(@(x) -2*H1*c1*(x-x1)*exp(-c1*(x - x1).^2) + 2*H2*c2*(x-x2)*exp(-c2*(x - x2).^2), [x0 Xf]);

%u   = chebfun(@(x) -eta(x)/sqrt(x), [0 L]);
%u_prime = chebfun(@(x) -( ( eta(x)/(2*sqrt(x)) + sqrt(x)*eta_prime(x) )/x ), [0 L]);  %needs to change

u_0 = @(x) 0.0;
u_prime = chebfun(@(x) 0.00, [x0 Xf]);

%Nicolsky 2018 solution
[ana_eta, ana_u] = catalina_transform(true, 1000); %data projection

%comparison
%[eta_diff, eta_l2] = stat(num_eta, ana_eta, 1000);
%[u_diff, u_l2] = stat(num_u, ana_u, 1000);
%[eta1_diff, eta1_l2] = stat(num_eta, cat_eta, 1000);
%[u1_diff, u1_l2] = stat(num_u, cat_u, 1000);

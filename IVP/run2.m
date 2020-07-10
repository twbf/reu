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
global eta_0 eta_prime u_0 u_prime td t0 Tf x0 Xf upS lowS upL lowL K

%%% bathymetry and time

td = 10.0/10.0; %slope of bathymetry

t0 = 0; %initial time
Tf = 3; %final time (10 in Catalina 1 experiement)

x0 = 0;
Xf = 7;

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

u_0   = chebfun(@(x) -2.*(sqrt(x+eta_0(x)+3)-sqrt(x+3)), [0 30]);
u_prime = diff(u_0);


%u_prime = chebfun(@(x) -( ( eta(x)/(2*sqrt(x)) + sqrt(x)*eta_prime(x) )/x ), [0 L]);  %needs to change

%u_0 = @(x) 0.1*exp(-(x-3).^2);
%u_prime = chebfun(@(x) -0.01*2*(x-3)*exp(-(x-3).^2), [x0 Xf]);

%Nicolsky 2018 solution

% integration parameters
lowS = x0;  %lower bound on s
upS = xF;   %upper bound on s
lowL = t0;  %lower bound on lambda
upL = tF;   %upper bound on lambda

K = 20; % upper bound on integration parameter

[ana_eta, ana_u] = CG_transform(true, 1000, true, false);

 t0 = 0; %initial time
 Tf = 0.9; %final time (10 in Catalina 1 experiement)

 x0 = 1;
 Xf = 6;
%numeric solution
[num_eta, num_u] = run_num();


u_0   = @(x) 0;

[num1_eta, num1_u] = run_num();


%catalina 1 solution
%[cat_eta, cat_u] = run_num();

disp(' ');
disp('Comparing');

%comparison
[diff, l2] = stat(num_eta, ana_eta, 1000);
%[u_diff, u_l2] = stat(num_u, ana_u, 1000);
%[eta1_diff, eta1_l2] = stat(num_eta, cat_eta, 1000);
%[u1_diff, u1_l2] = stat(num_u, cat_u, 1000);

disp(' ');
disp('Saving and Ploting');

%display
save('test1');
Plot_xt(num_eta, 500, 1, 'Eta Numerical Solution (Deny FV)', 'x', 't', 'eta');
Plot_xt(num1_eta, 500, 2, 'Eta Anaylytic', 'x', 't', 'eta');
% Plot_xt(eta_diff, 500, 3, 'Eta Diff', 'x', 't', 'diff')
%
%Plot_xt(num_u, 500, 3, 'U Numerical Solution (Deny FV)', 'x', 't', 'u')
Plot_xt(ana_eta, 500, 3, 'U Anaylytic', 'x', 't', 'u')
Plot_xt(diff, 500, 4, 'U Diff', 'x', 't', 'diff')

%Plot_xt(cat_u, 500, 7, 'U Catalina 1', 'x', 't', 'eta')
%Plot_xt(cat_eta, 500, 8, 'Eta Catalina1', 'x', 't', 'u')
%Plot_xt(eta1_diff, 500, 9, 'Eta1 Diff', 'x', 't', 'diff')

figure(5)
plot(linspace(t0, Tf, 1000), l2);

%figure(11)
%plot(linspace(t0, Tf, 1000), u_l2);

%figure(12)
%plot(linspace(t0, Tf, 1000), eta1_l2);

%figure(13)
%plot(linspace(t0, Tf, 1000), u1_l2);

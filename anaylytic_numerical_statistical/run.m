
%notes: look at all the ranges on the chgebfuns

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

td = 1.0/10.0; %slope of bathymetry

t0 = 0; %initial time
Tf = 10; %final time (10 in Catalina 1 experiemen)

x0 = 0;
Xf = 10;

%%%initial conditions

H1 = 0.006;
H2 = 0.018;
c1 = 0.4444;
c2 = 4.0;
x1 = 4.1209;
x2 = 1.6384;

eta_0 = chebfun(@(x) H1*exp(-c1*(x - x1).^2) - H2*exp(-c2*(x - x2).^2), [x0 Xf]);
eta_prime = chebfun(@(x) -2*H1*c1*(x-x1)*exp(-c1*(x - x1).^2) + -2*H2*c2*(x-x2)*exp(-c2*(x - x2).^2), [x0 Xf]);

%u   = chebfun(@(x) -eta(x)/sqrt(x), [0 L]);
%u_prime = chebfun(@(x) -( ( eta(x)/(2*sqrt(x)) + sqrt(x)*eta_prime(x) )/x ), [0 L]);  %needs to change

u_0 = chebfun(@(x) 0, [x0 Xf]);
u_prime = chebfun(@(x) 0, [x0 Xf]);


%numeric solution

[num_eta, num_u] = run_num();

%Nicolsky 2018 soultion

%catalina 1 solution

%comparison

%display

num_x = 100;
num_t = 100;

t_list = linspace(t0, Tf, num_t);
x_list = linspace(x0, Xf, num_x);

xx = zeros(num_x, num_t);
tt = zeros(num_x, num_t);
hh = zeros(num_x, num_t);

for t=1:num_t
  for x=1:num_x
    xx(x,t) = t_list(t);
    tt(x,t) = x_list(x);
    hh(x,t) = num_eta(x_list(x),t_list(t));
  end
end



figure(1)
mesh(tt,xx,hh)
title(['Numerical Solution (Deny FV)']);
xlabel('x');
ylabel('t');

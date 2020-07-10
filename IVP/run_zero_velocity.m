
clear
close all
format longE

%%%libraries we use
addpath('numerical/');
addpath('catalina1/');
addpath('anaylytic/');

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


%%% catalina 1 solution
%[cat_eta, cat_u] = catalina_transform(true, 1000);

%%% numeric solution
[num_eta, num_u] = run_num();

%%% Nicolsky 2018 solution
% integration parameters
lowS = x0;  %lower bound on s
upS = xF;   %upper bound on s
lowL = t0;  %lower bound on lambda
upL = tF;   %upper bound on lambda

K = 20; % upper bound on integration parameter

[ana_eta, ana_u] = CG_transform(false, 4000, false, true);

%comparison
[eta_diff, eta_l2] = stat(num_eta, ana_eta, 1000);
[u_diff, u_l2] = stat(num_u, ana_u, 1000);
[eta1_diff, eta1_l2] = stat(num_eta, cat_eta, 1000);
[u1_diff, u1_l2] = stat(num_u, cat_u, 1000);

%display
Plot_xt(num_eta, 500, 1, 'Eta Numerical Solution (Deny FV)', 'x', 't', 'eta')
Plot_xt(ana_eta, 500, 2, 'Eta Anaylytic', 'x', 't', 'u')
Plot_xt(eta_diff, 500, 3, 'Eta Diff', 'x', 't', 'diff')

Plot_xt(num_u, 500, 4, 'U Numerical Solution (Deny FV)', 'x', 't', 'eta')
Plot_xt(ana_u, 500, 5, 'U Anaylytic', 'x', 't', 'u')
Plot_xt(u_diff, 500, 6, 'U Diff', 'x', 't', 'diff')

Plot_xt(cat_u, 500, 7, 'U Catalina 1', 'x', 't', 'eta')
Plot_xt(cat_eta, 500, 8, 'Eta Catalina1', 'x', 't', 'u')
Plot_xt(eta1_diff, 500, 9, 'Eta1 Diff', 'x', 't', 'diff')

figure(10)
plot(linspace(t0, Tf, 1000), eta_l2);

figure(11)
plot(linspace(t0, Tf, 1000), u_l2);

figure(12)
plot(linspace(t0, Tf, 1000), eta1_l2);

figure(13)
plot(linspace(t0, Tf, 1000), u1_l2);

% run
% j0...
% j1...
% cos...
% sin...
% 2...
% 1...
% 3...
% Error using chebfun3/subsref (line 46)
% Unrecognized inputs.
%
% Error in catalina1>@(sigma,lambda)sum(xi^2*Phi(sigma)*inner_e_sl_1(lambda,sigma,xi)) (line 51)
%     eta_sl_1 = chebfun2( @(sigma, lambda) sum( xi^2*Phi(sigma) * inner_e_sl_1(lambda, sigma, xi)), [0.01 .1 0
%     .1], 'vectorize');
%
% Error in chebfun2/constructor>evaluate (line 409)
%             vals(jj, kk) = op( xx( 1, kk) , yy( jj, 1 ) );
%
% Error in chebfun2/constructor (line 76)
%     vals = evaluate(op, xx, yy, vectorize);
%
% Error in chebfun2 (line 82)
%             f = constructor(f, varargin{:});
%
% Error in catalina1 (line 51)
%     eta_sl_1 = chebfun2( @(sigma, lambda) sum( xi^2*Phi(sigma) * inner_e_sl_1(lambda, sigma, xi)), [0.01 .1 0
%     .1], 'vectorize');
%
% Error in run (line 46)
% [cat_eta, cat_u] = catalina1();

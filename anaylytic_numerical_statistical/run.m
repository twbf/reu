
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

td = 0.4/10.0; %slope of bathymetry

t0 = 0; %initial time
Tf = 10; %final time (10 in Catalina 1 experiement)

x0 = -2;
Xf = 10;

%%%initial conditions

H1 = 0.006;
H2 = 0.018;
c1 = 0.4444;
c2 = 4.0;
x1 = 4.1209;
x2 = 1.6384;

eta_0 = @(x) H1*exp(-c1*(x - x1).^2) - H2*exp(-c2*(x - x2).^2);
eta_prime = chebfun(@(x) -2*H1*c1*(x-x1)*exp(-c1*(x - x1).^2) + -2*H2*c2*(x-x2)*exp(-c2*(x - x2).^2), [x0 Xf]);

%u   = chebfun(@(x) -eta(x)/sqrt(x), [0 L]);
%u_prime = chebfun(@(x) -( ( eta(x)/(2*sqrt(x)) + sqrt(x)*eta_prime(x) )/x ), [0 L]);  %needs to change

%u_0 = @(x) 0.00001*exp(x-5);
u_0 = @(x) 0;
u_prime = chebfun(@(x) 0, [x0 Xf]);

%catalina 1 solution
%[cat_eta, cat_u] = catalina1();

%numeric solution
[num_eta, num_u] = run_num();
td = 1.0/10.0;
[num_eta_0, num_u_0] = run_num();

%Nicolsky 2018 solution
%[ana_eta, ana_u] = CG_transform(false, 4000, false, true);

%comparison
[an_diff, an_l2] = stat(num_eta, num_eta_0, 1000);

%display
Plot_xt(num_eta, 500, 1, 'Eta Numerical Solution (Deny FV)', 'x', 't', 'eta')
Plot_xt(num_eta_0, 500, 2, 'U Numerical Solution (Deny FV)', 'x', 't', 'u')
Plot_xt(an_diff, 500, 3, 'Diff Numerical Solution (Deny FV)', 'x', 't', 'diff')

figure(4)
plot(linspace(t0, Tf, 1000), an_l2);

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

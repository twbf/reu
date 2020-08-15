%---------------RUNSCRIPT--------------------------%

addpath('IV_FVM/')
addpath('data_projection/')
addpath('FD_Hodograph/')

%%% Global variables:
global t0 Tf x0 Xf
global t_res x_res
global eta_0 u_0 g td

% Bottom slope:
td = 10.0/10.0;

% initial time and final time (domain of t)
t0 = 0;
Tf = 1;

% domain of x
x0 = -2;
Xf = 10;

% resolution parameters
t_res = 100;
x_res = 2000;

%% Initial condition parameters:
H1 = 0.006;
H2 = 0.018;
c1 = 0.4444;
c2 = 4.0;
x1 = 4.1209;
x2 = 1.6384;

eta_0 = @(x) H1*exp(-c1*(x - x1).^2) - H2*exp(-c2*(x - x2).^2);
u_0   = @(x) 0;

[eta1, u1] = BCpull();   % solution via FVM
Psi = order_n_BC_proj(eta1, u1, 1);   % Data Projection

figure(3);
plot(Psi(1,:));

figure(4);
plot(Psi(2,:));

[phi, psi] = HodoSolve(Psi); % hodograph Solver

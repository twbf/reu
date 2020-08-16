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
Tf = 1.7;

% domain of x
x0 = -2;
Xf = 10;

% resolution parameters
t_res = 5000;
x_res = 6000;

%% Initial condition parameters:
H1 = 0.006;
H2 = 0.018;
c1 = 0.4444;
c2 = 4.0;
x1 = 4.1209+2;
x2 = 1.6384+2;

g = 9.81;

eta_0 = @(x) H1*exp(-c1*(x - x1).^2) - H2*exp(-c2*(x - x2).^2);
u_0   = @(x) 0;

[eta1, u1, eta_fvm, u_fvm] = BCpull();   % solution via FVM
Psi = order_n_BC_proj(eta1, u1, 4);   % Data Projection

  figure(3);
  plot(Psi(1,:));

  figure(4);
  plot(Psi(2,:));

test = [eta1; u1/sqrt(g)];

[phi, psi] = HodoSolve(Psi); % hodograph Solver

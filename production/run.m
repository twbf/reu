
%   GENERAL INITIAL CONDITION ANALYTIC VS NUMERICAL SOLUTION COMPARISON
%
% Works with a given IC on a given domain. Produces a comparison of the analytic
% solution
%
% Notes:
%    a(k) and b(k) in "Nicolsky_2018/fast_hankel" must converge to zero
%       before the end of the k domain. Please adjust both the resolution
%       and domain for this to effectively work
%
%

clear
close all
format longE

addpath('data_projection/');
addpath('IV_FVM/');
addpath('Nicolsky_2018/');

global eta_0 eta_prime u_0 u_prime td g t0 Tf x0 Xf %physical variables
global x_res t_res %resolution
global x k la s  %variables

disp('initializing variables ... ');

H1 = 0.006;
H2 = 0.018;
c1 = 0.4444;
c2 = 4.0;
x1 = 4.1209;
x2 = 1.6384;

% initial time and final time (domain of t)
t0 = 0;
Tf = 5;

% doman of x
x0 = -2;
Xf = 10;

% defining IC - chebfuns are used such that eta_prime and u_prime
% can be defeined simply for any function
eta_0 = chebfun(@(x) H1*exp(-c1*(x - x1).^2) - H2*exp(-c2*(x - x2).^2), [x0, Xf]);
%eta_0 = chebfun(@(x) 0, [x0, Xf]);
eta_prime = diff(eta_0);

u_0   = chebfun(@(x) -0.03*sin(3*x)*exp(-0.5*(x-5)^2), [x0 Xf]);
%u_0   = chebfun(@(x) 0, [x0 Xf]);
u_prime = diff(u_0);

td = 10.0/10.0; % slope
g = 9.81; % gravity acceleration


x_res = 2000; % number of points in the x domain - also used for s and k
t_res = 100; % number of points in the t domain - also used for lamda


%setting variables
x = linspace(x0,Xf,x_res);
t = linspace(t0,Tf,t_res);
k = linspace(0,70,x_res);
la = linspace(t0,Tf*sqrt(g),t_res);
s = linspace(0,Xf,x_res);


%computing solutions
[eta_analytic u_analytic] = fast_hankel(3);
[eta_fvm u_fvm] = run_num();


%post processing

disp('post processing eta.....');

num = zeros(t_res, x_res);
ana = zeros(t_res, x_res);

for i = 1:t_res

  ana(i,:) = eta_analytic(x, repmat(t(i), 1, x_res));
  num(i,:) = eta_fvm(:, i)';

  %making it zero on the other side of the shore
  max = 0;
  for j = 1:x_res
    if num(i, j) + td*x(j) < 0
      num(i, j) = NaN;
      max = j;
    end
    if ana(i,j) + td*x(j) < 0
      ana(i,j) = NaN;
      max = j;
    end
  end
  stat_norm(i) = norm(ana(i, max+1:end) - num(i, max+1:end));
end

disp('displaying .... ')

figure(5);
%disp(max)
for i = 1:t_res
  plot(x, num(i,:) ), hold on;
  plot(x, ana(i,:) ), hold on;
  plot(x, -td*x), hold off;
  axis([x0 Xf -0.05 0.05])
  pause(0.01);
end

figure(5);
subplot(4,1,1);
plot(x, num(1,:)), hold on;
% plot(x, -td*x), hold on;
plot(x, ana(1,:)), hold off;

subplot(4,1,2);
plot(x, num(20,:)), hold on;
%plot(x, -td*x), hold on;
plot(x, ana(20,:)), hold off;

subplot(4,1,3);
plot(x, num(60,:)), hold on;
%plot(x, -td*x), hold on;
plot(x, ana(60,:)), hold off;

subplot(4,1,4);
plot(x, stat_norm);




figure(6);
plot(x, stat_norm);

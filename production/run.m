clear
close all
format longE

addpath('data_projection/');
addpath('IV_FVM/');
addpath('Nicolsky_2018/');

global eta_0 eta_prime u_0 u_prime td g t0 Tf x0 Xf
global x_res t_res %resolution
global x k la s  %variables

disp('initializing variables ... ');

H1 = 0.006;
H2 = 0.018;
c1 = 0.4444;
c2 = 4.0;
x1 = 4.1209;
x2 = 1.6384;

t0 = 0;
Tf = 5;

x0 = -2;
Xf = 10;

%eta_0 = chebfun(@(x) H1*exp(-c1*(x - x1).^2) - H2*exp(-c2*(x - x2).^2), [x0, Xf]);
eta_0 = chebfun(@(x) 0.01*sin(10*x)*exp(-0.3*(x-5)^2), [x0, Xf]);
eta_prime = diff(eta_0);

%u_0   = chebfun(@(x) -0.03*sin(3*x)*exp(-0.5*(x-5)^2), [x0 Xf]);
u_0   = chebfun(@(x) 0, [x0 Xf]);
u_prime = diff(u_0);

td = 10.0/10.0; %slope
g = 9.81; %gravity acceleration


x_res = 2000;
t_res = 500;

x = linspace(x0,Xf,x_res);
t = linspace(t0,Tf,t_res);
k = linspace(0,70,x_res);
la = linspace(t0,Tf*sqrt(g),t_res);
s = linspace(0,Xf,x_res);

%running solutions

disp('analytic solution... ');
[eta_analytic u_analytic] = fast_hankel();

[eta_fvm u_fvm] = run_num();


%post processing
close all
figure(1);
x = linspace(x0,Xf,x_res);
for i = 1:t_res

  ana = eta_analytic(x, repmat(t(i), 1, x_res));
  num = eta_fvm(:, i)';

  %making it zero on the other side of the shore
  max = 0;
  for j = 1:x_res
    if num(j) + td*x(j) < 0
      num(j) = NaN;
      max = j;
    end
    if ana(j) + td*x(j) < 0
      ana(j) = NaN;
      max = j;
    end
  end
  %disp(max)
  stat_norm(i) = norm(ana(max+1:end) - num(max+1:end));
  plot(x, num ), hold on;
  plot(x, ana ), hold on;
  plot(x, -td*x), hold off;
  axis([x0 Xf -0.05 0.05])
  pause(0.01);
end

figure(7)
plot(stat_norm)

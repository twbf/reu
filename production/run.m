clear
close all
format longE

addpath('data_projection/');
addpath('IV_FVM/');
addpath('Nicolsky_2018/');

global eta_0 eta_prime u_0 u_prime td t0 Tf x0 Xf
global x_res t_res %resolution
global x k la s  %variables

H1 = 0.006;
H2 = 0.018;
c1 = 0.4444;
c2 = 4.0;
x1 = 4.1209;
x2 = 1.6384;

t0 = 0;
Tf = 3;

x0 = -2;
Xf = 10;

eta_0 = chebfun(@(x) H1*exp(-c1*(x - x1).^2) - H2*exp(-c2*(x - x2).^2), [x0, Xf]);
eta_prime = diff(eta_0);

u_0   = chebfun(@(x) 0, [x0 Xf]);
u_prime = diff(u_0);

td = 10.0/10.0;

x_res = 3000;
t_res = 800;

x = linspace(x0,Xf,x_res);
t = linspace(t0,Tf,t_res);
k = linspace(0,30,x_res);
la = linspace(t0,Tf*sqrt(9.81),t_res);
s = linspace(0,Xf,x_res);

%running solutions

[eta_analytic u_analytic] = fast_hankel();

[eta_fvm u_fvm] = run_num();


%post processing
figure(6);
x = linspace(x0,Xf,x_res);
for i = 1:t_res

  ana = eta_analytic(x, repmat(t(i), 1, x_res));
  num = eta_fvm(:, i)';

  %making it zero on the other side of the shore
  for j = 1:x_res
    if num(j) + x(j) < 0
      num(j) = 0;
    end
    if ana(j) + x(j) < 0
      ana(j) = 0;
    end
  end

  stat_norm(i) = norm(ana - num);
  plot(x, num ), hold on;
  plot(x, ana ), hold on;
  plot(x, -x), hold off;
  axis([x0 Xf -0.03 0.03])
  pause(0.001);
  disp(i)
end

figure(7)
plot(stat_norm)

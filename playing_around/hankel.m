clear
close all
format longE
% just testing out fast foureir transform

num_x = 1000;
num_k = 1000;
num_s = 1000;
num_la = 50;

x = linspace(0,10,num_x);
k = linspace(0,30,num_k);
la_list = linspace(0,2.5,num_la);
s = linspace(0,10,num_s);


g = zeros(size(x));
f = zeros(size(x));


%for i=1:num_x
%  s_s = ss(x(i));
%  p = proj(x(i));
%  g(i) = [0 1]*p;
  %f(i) = [1 0]*p*s_s.^(1/2);
%end


%test1 = 2*k.*ihat(sss(x, true), sqrt(x), k, 0);
%test  = 2*k.*ifht(sss(x, true), x, k, 0);

%disp(test)
%figure(1)
%plot(test);



b = -2*k.*ihat(sss(x, false), sqrt(x), 2*k, 1);

%a = 20*k.*ihat(sss(x, true), sqrt(x), 2*k, 0);

load('zero_initial_k30', 'a');

a = a(k);

figure(1);
plot(k, a);

figure(2);
plot(k, b);

psi = @(la) a.*cos(la*k) + b.*sin(la*k);
phi = @(la) a.*sin(la*k) - b.*cos(la*k);




%to display the animation
%for i=1:num_la
%  disp(i)
%  plot(x, hat(psi(la_list(i)),k));
%  pause(0.05);
%end

Psi = zeros(num_la, num_k);
Phi = zeros(num_la, num_k);

for i=1:num_la
  disp(i)
  Psi(i, :) = hat(psi(la_list(i)), k, sqrt(s), 0);
  Phi(i, :) = s.^(-1/2).*hat(phi(la_list(i)),  k, sqrt(s), 1); %replaced h with the correct phi
end


figure(3);
mesh(Psi);

figure(4);
mesh(Phi);

eta = zeros(num_la, num_x);
u = zeros(num_la, num_x);
xx = zeros(num_la, num_x);
tt = zeros(num_la, num_x);

x = s;

for i=1:num_x
  %CG transform
  u(:,i) = Phi(:,i);
  eta(:,i) = Psi(:,i) - u(:,i).^2/2;
  xx(:,i) = x(i)*repmat(1,num_la,1) - eta(:,i);
  tt(:,i) = u(:,i) + la_list';

  %deminsionalizing
  u(:,i) = u(:,i)*sqrt(9.81);
  tt(:,i) = tt(:,i)./sqrt(9.81);
end

figure(5);
mesh(tt,xx,eta);

figure(6);
scatter(tt,xx);

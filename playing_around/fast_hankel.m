clear
close all
format longE
% just testing out fast foureir transform

num_x = 1000;
num_k = 1000;
num_la = 100;

x = linspace(0,30,num_x);
k = linspace(0,70,num_k);
la_list = linspace(0,3,num_la);
eta_0 = @(x) 0.1*exp(-(x-5).^2);
eta_prime = @(x) -2*0.1*(x-5).*exp(-(x-5).^2);

H1 = 0.006;
H2 = 0.018;
c1 = 0.4444;
c2 = 4.0;
x1 = 4.1209 +4;
x2 = 1.6384 +4;

%eta_0 = chebfun(@(x) H1*exp(-c1*(x - x1).^2) - H2*exp(-c2*(x - x2).^2), [-3, 50]);
%eta_prime = diff(eta_0);



%u_0   = chebfun(@(x) -eta_0(x)/sqrt(x+eps), [0 30]);
%u_prime = diff(u_0);
u_0 = @(x) 0;
u_prime = @(x) 0;

s = @(x) x + eta_0(x);
A = @(x) [0 1; s(x) 0];
B = [0 0; 1 0];

D = @(x) (1+eta_prime(x))*eye(2) + u_prime(x)*A(x);

phi0 = @(x) [u_0(x) ; eta_0(x)+(u_0(x).^2)/2];
phi0_prime = @(x) [u_prime(x); eta_prime(x)+u_0(x)*u_prime(x)];

proj = @(x) phi0(x) + u_0(x).*(u_prime(x).*A(x)*inv(D(x))*B*phi0(x) - B*phi0(x) - A(x)*inv(D(x))*phi0_prime(x));

g = zeros(size(x));
f = zeros(size(x));

for i=1:num_x
  p = proj(x(i));
  g(i) = [0 1]*p;
  f(i) = [1 0]*p*x(i).^(1/2);
end

ss = @(x) cos(x);
[test t1 t2 t3 t4] = fht(g, 10, 100, [1,2], 7, 5);

disp(test(2));

plot(t4)

disp(test(1));

b = -2*k.*hat(f,k, pi/numel(f)*(0:numel(f)-1),1 );
a = 2*k.*hat(g,k);

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
  Psi(i, :) = hat2(psi(la_list(i)),k)./(num_k.^2);
  Phi(i, :) = x.^(-1/2).*hat2(phi(la_list(i)), k, pi/numel(f)*(0:numel(f)-1),1)./(num_k.^2);
end

figure(3);
mesh(Psi);

figure(4);
mesh(Phi);

eta = zeros(num_la, num_x);
u = zeros(num_la, num_x);
xx = zeros(num_la, num_x);
tt = zeros(num_la, num_x);

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

clear
close all
format longE
% just testing out fast foureir transform

num_x = 2000;
num_k = 2000;
num_s = 2000;
num_la = 2;

x = linspace(0.1,80,num_x);
k = linspace(0.1,100,num_k);
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


%disp(test)
%figure(1)
%plot(test);



b = -2*k.*ihat(sss(x, false), sqrt(x), 2*k, 1);

%a_fht = ifht(sss(x, true), k , sqrt(x),  0);

%tda = @(vx) interp1(x, sss(x, true), vx);
a_fht  = ifht(sss(x, true), k, x);




%a_ihat = 2*k.*ihat(sss(x, true), sqrt(x), 2*k, 0);

figure(1);
plot(k, a_fht)


stop

load('zero_initial_k30', 'a');

upS = 10;
K = 30;
upL = 10;

y = chebfun('x', [0 K]);

disp('j0...')
j0 = chebfun(@(kx) besselj(0.0, kx), [0.0 max( [2.0*30*sqrt(upS)])]);

disp('cos...')
Cos = chebfun(@(lk) cos(lk), [0 upL*K], 'vectorize');

disp('sin...')
Sin = chebfun(@(lk) sin(lk), [0 upL*K], 'vectorize');

psi = chebfun(@(s) sum( ( a(y)*Cos(1*y)) * j0(2.0*y*sqrt(s)) ), [0 10], 'vectorize');
%psi = chebfun2(@(s,la) sum( ( a(k)*Cos(la*k) + b(k)*Sin(la*k) ) * j0(2.0*k*sqrt(s)) ), [0 upS 0 upL], 'vectorize');

%psi_fht = @(k) a(k).*cos(1*k);


a = a(k);

psi_fht = @(vk)  interp1(k, a , vk).*cos(1*vk);


figure(2);
plot(a_ifht);




psi = @(la) a.*cos(la*k) + b.*sin(la*k);
psi1 = @(la) a_ihat.*cos(la*k) + b.*sin(la*k);
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
  Psi(i, :) = hat(psi(la_list(i)), k, 2*sqrt(s), 0);
  Psi_1(i, :) = hat(psi1(la_list(i)), k, 2*sqrt(s), 0);
  [Psi_fht(i, :) t3]= fht(psi_fht, 30, 7, 0, 100, 99);
  %Psi_fht(i, :) = fht(a, k, k, k);
  Phi(i, :) = s.^(-1/2).*hat(phi(la_list(i)),  k, 2*sqrt(s), 1); %replaced h with the correct phi
end


figure(3);
mesh(Psi);

%figure(4);
%mesh(Phi);

figure(4);
plot(s, Psi(1, :))

figure(5);
plot(t3.^2./4, Psi_fht(1, :))

eta = zeros(num_la, num_x);
u = zeros(num_la, num_x);
xx = zeros(num_la, num_x);
tt = zeros(num_la, num_x);

x = s;

for i=1:num_x
  %CG transform
  u(:,i) = Phi(:,i);
  eta(:,i) = Psi(:,i) - u(:,i).^2/2/40;
  xx(:,i) = x(i)*repmat(1,num_la,1) - eta(:,i);
  tt(:,i) = u(:,i) + la_list';

  %deminsionalizing
  u(:,i) = u(:,i)*sqrt(9.81);
  tt(:,i) = tt(:,i)./sqrt(9.81);
end

figure(6);
mesh(tt,xx,eta);

figure(7);
scatter(tt,xx);

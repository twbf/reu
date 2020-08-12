clear
close all
format longE
% just testing out fast foureir transform

num_x = 1000;
num_k = 1000;
num_s = 1000;
num_la = 80;

x_end = 10;

x = linspace(0,x_end,num_x);
k = linspace(0,30,num_k);
la_list = linspace(0,10,num_la);
s = linspace(0,10,num_s);

x_density = num_x/x_end;

b = -2*k.*ihat(sss(x, false), sqrt(x), 2*k, 1)./x_density;
a = 2*k.*ihat(sss(x, true), sqrt(x), 2*k, 0)./x_density;

figure(1);
plot(k, a);

%load('zero_initial_k30', 'a');
%a = a(k);

figure(2);
plot(k, b);


for i=1:num_la
  disp(i)

  psi = @(vk)  interp1(k, a , vk).*cos(la_list(i)*vk) + interp1(k, b , vk).*sin(la_list(i)*vk);
  phi = @(vk)  interp1(k, a , vk).*sin(la_list(i)*vk) - interp1(k, b , vk).*cos(la_list(i)*vk);

  [Psi(i, :) r_psi] = fht(psi, 30, 7, 0, 20, 15);
  [Phi(i,:) r_phi] = fht(phi, 30, 7, 1, 20, 15); % needsa to include s^(-1/2)

  %need to deal with scalling
  Psi(i, :) = Psi(i, :)./(2*pi);
  Phi(i, :) = r_phi.^(1/2).*Phi(i, :)./(2*pi);

end

r_size = size(Psi); %same for phi and psi

figure(3);
mesh(repmat(r_psi.^2./4, num_la, 1), repmat(la_list.', 1, r_size(2)), Psi);

figure(4);
mesh(repmat(r_phi.^2./4, num_la, 1), repmat(la_list.', 1, r_size(2)), Phi);


eta = zeros(num_la, r_size(2));
u = zeros(num_la, r_size(2));
xx = zeros(num_la, r_size(2));
tt = zeros(num_la, r_size(2));

x = s;

for i=1:r_size(2)
  %CG transform
  u(:,i) = Phi(:,i);
  eta(:,i) = Psi(:,i) - u(:,i).^2/2;
  xx(:,i) = r_phi(i).^2./4 - eta(:,i);
  tt(:,i) = u(:,i) + la_list';

  %deminsionalizing
  %u(:,i) = u(:,i)*sqrt(9.81);
  %tt(:,i) = tt(:,i)./sqrt(9.81);
end

figure(5);
mesh(tt,xx,eta);

load('Denys_FVM/fvm_test', 'eta_fvm');

s_tt = reshape(tt, [num_la*r_size(2), 1]);
s_xx = reshape(xx, [num_la*r_size(2), 1]);
s_eta = reshape(eta, [num_la*r_size(2), 1]);
test = scatteredInterpolant(s_xx, s_tt,s_eta);

figure(6);
x_plot = linspace(-1, 10, 8000);
t_plot = linspace(0, 10, 800);
pause(4);
for i = 1:800

  ana = test(x_plot, repmat(t_plot(i), 1, 8000));
  num = eta_fvm(:, i)';

  disp(size(eta_fvm));

  for j = 1:8000
    if num(j) == 0
      ana(j) = 0;
    end
  end

  stat_norm(i) = norm(ana - num);

  %plot(x_plot, ana ), hold on;
  %%plot(x_plot, -x_plot), hold off;
  %axis([-1 10 -0.03 0.03])
  %pause(0.001);
  disp(i)
end

figure(7)
plot(stat_norm)

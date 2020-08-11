clear
close all
format longE
% just testing out fast foureir transform

num_x = 8000;
num_k = 8000;
num_s = 8000;
num_la = 800;

x = linspace(0.1,10,num_x);
k = linspace(0.1,30,num_k);
la_list = linspace(0,10,num_la);
s = linspace(0,10,num_s);

%b = -2*k.*ihat(sss(x, false), sqrt(x), 2*k, 1);
a1 = k.*ifht(sss(x, true), x, k, 0);
%a = fht(sss(x, true), 30, 7, 0);

figure(1);
plot(a1);

load('zero_initial_k30', 'a');
a = a(k);

figure(2);
plot(k, a-a1);

stop



for i=1:num_la
  disp(i)

  psi = @(vk)  interp1(k, a , vk).*cos(la_list(i)*vk) + interp1(k, b , vk).*sin(la_list(i)*vk);
  phi = @(vk)  interp1(k, a , vk).*sin(la_list(i)*vk) - interp1(k, b , vk).*cos(la_list(i)*vk);

  [Psi(i, :) r_psi] = fht(psi, 30, 7, 0, 20, 15);
  [Phi(i,:) r_phi] = fht(phi, 30, 7, 1, 20, 15); % needsa to include s^(-1/2)

  %need to deal with scalling
  Psi(i, :) = Psi(i, :)./6;
  Phi(i, :) = Phi(i, :)./6;
  %Phi(i, :) = s.^(-1/2).*hat(phi(la_list(i)),  k, 2*sqrt(s), 1); %replaced h with the correct phi
end

r_size = size(Psi) %same for phi and psi

figure(3);
mesh(repmat(r_psi.^2./4, num_la, 1), repmat(la_list.', 1, r_size(2)), Psi);
%for i = 1:num_la
%  plot(r_psi, Psi(i, :));
%  pause(0.05);
%end

figure(4);
mesh(repmat(r_phi.^2./4, num_la, 1), repmat(la_list.', 1, r_size(2)), Phi);
%for i = 1:num_la
%  plot(r_phi, Phi(i, :));
%  pause(0.05);
%end


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
%plot(eta_fvm(:, 1));


s_tt = reshape(tt, [num_la*r_size(2), 1]);
s_xx = reshape(xx, [num_la*r_size(2), 1]);
s_eta = reshape(eta, [num_la*r_size(2), 1]);
test = scatteredInterpolant(s_xx, s_tt,s_eta);
figure(6);
test(1,6);



x_plot = linspace(-1, 10, 8000);
t_plot = linspace(0, 10, 800);
pause(4);
for i = 1:800
  plot(x_plot, test(x_plot, repmat(t_plot(i), 1, 8000))), hold on;
  plot(x_plot, eta_fvm(:, i)'), hold on;
  plot(x_plot, -x_plot), hold off;
  axis([-1 10 -0.03 0.03])
  pause(0.0001);
end

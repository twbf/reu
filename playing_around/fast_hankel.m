clear
close all
format longE
% just testing out fast foureir transform

num_x = 4000;
num_k = 4000;
num_s = 4000;
num_la = 200;

x = linspace(0,10,num_x);
k = linspace(0,35,num_k);
la_list = linspace(0,10,num_la);
s = linspace(0,10,num_s);

b = -2*k.*ihat(sss(x, false), sqrt(x), 2*k, 1)./50;
a = 2*k.*ihat(sss(x, true), sqrt(x), 2*k, 0)./50;

figure(1);
plot(k, a);

figure(2);
plot(k, b);

for i=1:num_la
  disp(i)

  psi = @(vk)  interp1(k, a , vk).*cos(la_list(i)*vk) + interp1(k, b , vk).*sin(la_list(i)*vk);
  phi = @(vk)  interp1(k, a , vk).*sin(la_list(i)*vk) - interp1(k, b , vk).*cos(la_list(i)*vk);

  [Psi(i, :) r_psi] = fht(psi, 30, 7, 0, 20, 15);
  [Phi(i,:) r_phi] = fht(phi, 30, 7, 1, 20, 15); % needsa to include s^(-1/2)

  %need to deal with scalling
  Psi(i, :) = Psi(i, :);
  Phi(i, :) = Phi(i, :);
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



x_plot = linspace(-1, 10, 1000);
t_plot = linspace(0, 10, 200);
pause(4);
for i = 1:200
  plot(x_plot, test(x_plot, repmat(t_plot(i), 1, 1000))), hold on;
  plot(x_plot, eta_fvm(:, i)'), hold on;
  plot(x_plot, -x_plot), hold off;
  axis([-1 10 -0.2 0.2])
  pause(0.05);
  stop
end

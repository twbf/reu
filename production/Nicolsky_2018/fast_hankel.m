
%     FAST HANKEL SOLUTION TO 1 SPACIAL DEMENTION SWE
%
% We use an inverse hankel transfrom (not fast) for the computation of
% a(k) and b(k). Then we use a fast hankel transform for the compuatation of
% phi and psi. Finally, we use the CG transform to go back to original variables
% in (x,t) then return a scattered interpolant of eta and u.
%
% Refernce = [Nickolsky et. al. 2018]

function [eta u] = fast_hankel()

  disp('analytic solution... ');

  global x k la s  %variables
  global x_res t_res Xf g td %resolution

  x_density = x_res/Xf;

  %data_projection

  disp('    data_projection onto \lambda = 0... ');
  proj = order1_dp(s);

  global eta_0 u_0 %for comparing with data_projection

  figure(1);
  plot(s, proj(1, :)), hold on;
  plot(s, u_0(s));


  figure(2);
  plot(s, proj(2, :)), hold on;
  plot(s, eta_0(s)+u_0(s).^2/2), hold on;
  plot(s, eta_0(s));


  disp('    inverse hankel transform to compute a(k) and b(k)... ');

  %for the analytic solution to work properly both a and b need to be very
  %close to zero at the end of the domain of k
  a = 2*k.*ihat(proj(2, :) , sqrt(s), 2*k, 0)./x_density;
  b = -2*k.*ihat(proj(1, :) , sqrt(s), 2*k, 1)./x_density;

  figure(3);
  plot(k, a);

  figure(4);
  plot(k, b);


  disp('    fast hankel transform to compute psi and phi... ');

  for i=1:t_res

    psi_freq = @(vk)  interp1(k, a , vk).*cos(la(i)*vk) + interp1(k, b , vk).*sin(la(i)*vk);
    phi_freq = @(vk)  interp1(k, a , vk).*sin(la(i)*vk) - interp1(k, b , vk).*cos(la(i)*vk);

    [psi(i, :) r_psi] = fht(psi_freq, 30, 7, 0, 20, 15);
    [phi(i,:) r_phi] = fht(phi_freq, 30, 7, 1, 20, 15); % needs to include s^(-1/2)

    %need to deal with scalling
    psi(i, :) = psi(i, :)./(2*pi);
    phi(i, :) = r_phi.^(1/2).*phi(i, :)./(2*pi);
  end

  r_size = size(psi); %same for phi and psi


  %to display phi and psi
    %figure(3);
    %mesh(repmat(r_psi.^2./4, t_res, 1), repmat(la.', 1, r_size(2)), psi);

    %figure(4);
    %mesh(repmat(r_phi.^2./4, t_res, 1), repmat(la.', 1, r_size(2)), phi);

  disp('    backwards CG transform and demensionalization... ');

  eta = zeros(t_res, r_size(2));
  u = zeros(t_res, r_size(2));
  xx = zeros(t_res, r_size(2));
  tt = zeros(t_res, r_size(2));

  for i=1:r_size(2)
    %CG transform
    u(:,i) = phi(:,i);
    eta(:,i) = psi(:,i) - u(:,i).^2/2;
    xx(:,i) = r_phi(i).^2./4 - eta(:,i);
    tt(:,i) = u(:,i) + la';

    %deminsionalizing
    u(:, i) = u(:, i)*sqrt(g*td);
    eta(:, i) = eta(:, i);
    tt(:, i) = tt(:, i)/sqrt(td*g);
  end

  %to display eta and u
    %figure(5);
    %mesh(tt,xx,eta);

    %figure(5);
    %mesh(tt,xx,u);

  disp('    scattered interpolation of eta and u... ');

  s_tt = reshape(tt, [t_res*r_size(2), 1]);
  s_xx = reshape(xx, [t_res*r_size(2), 1]);
  s_eta = reshape(eta, [t_res*r_size(2), 1]);
  s_u = reshape(u, [t_res*r_size(2), 1]);

  eta = scatteredInterpolant(s_xx, s_tt, s_eta);
  u = scatteredInterpolant(s_xx, s_tt, s_u);

end

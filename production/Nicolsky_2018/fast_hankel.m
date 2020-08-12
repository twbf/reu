
function [eta u] = fast_hankel()

  global x k la s  %variables
  global x_res t_res Xf %resolution

  x_density = x_res/Xf;

  %data_projection
  proj = order1_dp(s);

  figure(1);
  plot(s, proj(1, :));

  figure(2);
  plot(s, proj(2, :));

  a = 2*k.*ihat(proj(2, :) , sqrt(s), 2*k, 0)./x_density;

  b = -2*k.*ihat(proj(1, :) , sqrt(s), 2*k, 1)./x_density;

  figure(3);
  plot(k, a);

  figure(4);
  plot(k, b);


  for i=1:t_res
    disp(i)

    psi_freq = @(vk)  interp1(k, a , vk).*cos(la(i)*vk) + interp1(k, b , vk).*sin(la(i)*vk);
    phi_freq = @(vk)  interp1(k, a , vk).*sin(la(i)*vk) - interp1(k, b , vk).*cos(la(i)*vk);

    [psi(i, :) r_psi] = fht(psi_freq, 30, 7, 0, 20, 15);
    [phi(i,:) r_phi] = fht(phi_freq, 30, 7, 1, 20, 15); % needs to include s^(-1/2)

    %need to deal with scalling
    psi(i, :) = psi(i, :)./(2*pi);
    phi(i, :) = r_phi.^(1/2).*phi(i, :)./(2*pi);
  end

  r_size = size(psi); %same for phi and psi

  figure(3);
  mesh(psi)

  figure(3);
  mesh(repmat(r_psi.^2./4, t_res, 1), repmat(la.', 1, r_size(2)), psi);

  figure(4);
  mesh(repmat(r_phi.^2./4, t_res, 1), repmat(la.', 1, r_size(2)), phi);

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
    u(:,i) = u(:,i)*sqrt(9.81);
    tt(:,i) = tt(:,i)./sqrt(9.81);
  end

  figure(5);
  mesh(tt,xx,eta);

  s_tt = reshape(tt, [t_res*r_size(2), 1]);
  s_xx = reshape(xx, [t_res*r_size(2), 1]);
  s_eta = reshape(eta, [t_res*r_size(2), 1]);
  s_u = reshape(u, [t_res*r_size(2), 1]);

  eta = scatteredInterpolant(s_xx, s_tt, s_eta);
  u = scatteredInterpolant(s_xx, s_tt, s_u);

end

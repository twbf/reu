
% Order 1 data_projection

% Refernce [Nicolsky et. al 2018]

function Phi = order1_dp(x)

  global eta_0 eta_prime u_0 u_prime

  ss = @(x) x + eta_0(x);
  A = @(x) [0 1; ss(x) 0];
  B = [0 0; 1 0];

  D = @(x) (1+eta_prime(x))*eye(2) + u_prime(x)*A(x);

  phi0 = @(x) [u_0(x) ; eta_0(x)+(u_0(x).^2)/2];
  phi0_prime = @(x) [u_prime(x); eta_prime(x)+u_0(x)*u_prime(x)];

  proj = @(x) phi0(x) + u_0(x).*(u_prime(x).*A(x)*inv(D(x))*B*phi0(x) - B*phi0(x) - A(x)*inv(D(x))*phi0_prime(x));


  x_num = size(x);

  Phi = zeros(2, x_num(2));

  for i = 1:x_num(2)
    Phi(:, i) = proj(x(i));
  end

end

function out = sss(inn, which)

  eta_0 = @(x) 0.1*exp(-(x-5).^2);
  eta_prime = @(x) -2*0.1*(x-5).*exp(-(x-5).^2);

  H1 = 0.006;
  H2 = 0.018;
  c1 = 0.4444;
  c2 = 4.0;
  x1 = 4.1209;
  x2 = 1.6384;


  eta_0 = chebfun(@(x) H1*exp(-c1*(x - x1).^2) - H2*exp(-c2*(x - x2).^2), [-2, 20]);
  eta_prime = diff(eta_0);

  %u_0   = chebfun(@(x) -eta_0(x)/sqrt(x+eps), [0 100]);
  %u_0   = chebfun(@(x) 0.01*sin(5*x), [0 100]);
  %u_prime = diff(u_0);
  u_0 = @(x) 0;
  u_prime = @(x) 0;

  ss = @(x) x + eta_0(x);
  A = @(x) [0 1; ss(x) 0];
  B = [0 0; 1 0];

  D = @(x) (1+eta_prime(x))*eye(2) + u_prime(x)*A(x);

  phi0 = @(x) [u_0(x) ; eta_0(x)+(u_0(x).^2)/2];
  phi0_prime = @(x) [u_prime(x); eta_prime(x)+u_0(x)*u_prime(x)];

  proj = @(x) phi0(x) + u_0(x).*(u_prime(x).*A(x)*inv(D(x))*B*phi0(x) - B*phi0(x) - A(x)*inv(D(x))*phi0_prime(x));

  out = zeros(size(inn));
  num = size(inn);
  for i = 1:num(2)
    p = proj(inn(i));

    if (which)
      out(i) = [0 1]*p;
    else
      out(i) = [1 0]*p*inn(i).^(1/2);
    end



    %g(i) = [0 1]*p;
    %f(i) = [1 0]*p*s_s.^(1/2);
  end

  %figure(10)
  %plot(out-eta_0(linspace(0,10,2000)))



end

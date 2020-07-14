function [eta_sl, u_sl] = catalina1()
    global  H1 H2 c1 c2 x1 x2 eta_0 eta_prime u_0 u_prime td t0 Tf x0 Xf

    % integration parameters
    upX = 50;                    % upper bound on xi
    upO = 50;                    % upper bound on omega
    upL = 3;                      % upper bound on lambda
    upS = 8;                     % upper bound on sigma

    %linear transform of x1, x2
    sigma1 = sqrt(x1*16);
    sigma2 = sqrt(x2*16);

    % eta(sigma, 0)
    % Equation (7) in catalina 1
    eta_sigma0 = chebfun(@(sigma) H1*exp(-(c1*(sigma^2-sigma1^2).^2)/256)...
    -H2*exp(-(c2*(sigma^2-sigma2^2).^2)/256), [0 upS]);

    % Plots initial wave function
    figure(1);
    plot(eta_sigma0);
    xlabel('$\sigma$','Interpreter','latex');
    ylabel('$\eta$','Interpreter','latex');

    % derivative of initial wave function
    eta_sigma0_prime = diff(eta_sigma0);

    % Plots the derivative of the initial wave function (respect to space)
    figure(3);
    plot(eta_sigma0_prime);
    xlabel('$\sigma$','Interpreter','latex');
    ylabel('$\eta^\prime$','Interpreter','latex');

   % Phi = chebfun(@(sigma) 4*eta_sigma0_prime(sigma)/sigma, [0 uS]);

   % Phi(sigma) given in equation (8)
   Phi = chebfun(@(sigma) -1/16*H1*c1*(sigma^2-sigma1^2)*exp(-1/256*c1*...
   (sigma^2-sigma1^2)^2) + 1/16*H2*c2*(sigma^2-sigma2^2)*exp(-1/256*c2*...
   (sigma^2-sigma2^2)^2) , [0 max(upS, upX)]);

   % Evaluates zeroth order bessel function for equation (4)
   disp('j0...')
   j0 = chebfun(@(kx) besselj(0.0, kx), [0 upO*upS]);
   % Evaluates firth order bessel function for equation (4)
   disp('j1...')
   j1 = chebfun(@(kx) besselj(1.0, kx), [0 max(upO*upX, upO*upS)]);
   % Evaluates cosine for equation (4)
   cosine = chebfun(@(lk) cos(lk), [0 upO*upL]);
   % Evaluates sine for equation (4)
   sine = chebfun(@(lk) sin(lk), [0 upO*upL]);
   % Evaluates xi for equation (4)
   xi = chebfun('x', [0 upX]);
   % Evaluates omega for equation (4)
   omega = chebfun('x', [0 upO]);
   % Evaluates inner integral for second term in equation (4)

   disp('u...')

   u_sl = chebfun2( @(sigma, lambda) sum(xi.^2*Phi(xi)*sum(j1(omega*sigma)./sigma*...
   j1(omega*xi)*sine(omega*lambda))), [0.001 upS 0 upL], 'vectorize');

   figure(2);
   plot(u_sl);
   xlabel('$\sigma$','Interpreter','latex');
   ylabel('$\lambda$','Interpreter','latex');
   zlabel('$u 2$','Interpreter','latex');

   %disp('u inner integral...')
   %in_int_term2 = chebfun2(@(sigma, lambda) sum(j1(omega*sigma)./sigma*...
   %j1(omega*xi)*sine(omega*lambda)), [0.001 upS 0 upL], 'vectorize');
   % Evaluates 'u' in equation (4)
   %u_sl = chebfun2( @(sigma, lambda) sum(xi.^2*Phi(xi)*in_int_term2...
   %(sigma, lambda)), [0.001 upS 0 upL], 'vectorize');
   %figure(3);
   %plot(u_sl);
   %xlabel('$\sigma$','Interpreter','latex');
   %ylabel('$\lambda$','Interpreter','latex');
   %zlabel('$u$','Interpreter','latex');


   % Evaluates inner integral for first term in equation (4)
   %disp('phi inner integral...')
   %in_int_term1 = chebfun2(@(sigma, lambda) sum(j0(omega*sigma)*j1(omega*xi)...
   %*cosine(omega*lambda)), [0.001 upS 0 upL], 'vectorize');
   % eta(sigma, lambda)

   disp('eta...')
   eta_sl = chebfun2(@(sigma,lambda) -1/4*sum(xi.^2*Phi(xi)*...
   sum(j0(omega*sigma)*j1(omega*xi)*cosine(omega*lambda)))-1/2*u_sl(sigma, lambda).^2,...
   [0.001 upS 0 upL], 'vectorize');

   figure(4);
   plot(eta_sl);
   xlabel('$\sigma$','Interpreter','latex');
   ylabel('$\lambda$','Interpreter','latex');
   zlabel('$\eta(\sigma,\lambda)$','Interpreter','latex');

   save('catalina1_phi_psi3');
end

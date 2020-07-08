%warning: look at the parameters to make sure all are good

function [eta_sl, u_sl] = catalina1(save)

    global H1 H2 c1 c2 x1 x2 eta_0 eta_prime u_0 u_prime td t0 Tf x0 Xf

    %integration Parameters
    Lx = 10; %upper bound on xi
    Lo = 10; %upper bound on omega
    Ll = 3; %upper bound on lambda
    Ls = 10; %upper bound on sigma

    %linear transform of x1, x2
    sigma1 = sqrt(x1*16);
    sigma2 = sqrt(x2*16);

    %eta(sigma, 0)
    eta_sigma_0 = chebfun(@(sigma) H1*exp( -(c1*(sigma^2 - sigma1^2).^2)/256 ) - H2*exp( -(c2*(sigma^2 - sigma2^2).^2)/256 ), [0 Ls]);

    figure(1)
    plot(eta_sigma_0)

    %eta_sigma_prime_0 = chebfun(@(sigma) diff(eta_sigma_0(sigma)), [0 Ls]);
    eta_sigma_prime_0 = chebfun(@(sigma) -4*H1*sigma*(c1*(sigma^2 - sigma1^2))/256*exp( -(c1*(sigma^2 - sigma1^2).^2)/256 ) + 4*H2*sigma*(c2*(sigma^2 - sigma2^2))/256*exp( -(c2*(sigma^2 - sigma2^2).^2)/256 ), [0 Ls]);

    figure(2)
    plot(eta_sigma_prime_0)

    %Phi = chebfun(@(sigma) 4*eta_sigma_prime_0(sigma)/sigma, [0 Ls]);
    Phi = chebfun(@(sigma) -1/16*H1*c1*(sigma^2-sigma1^2)*exp(-1/256*c1*(sigma^2-sigma1^2)^2) + 1/16*H2*c2*(sigma^2-sigma2^2)*exp(-1/256*c2*(sigma^2-sigma2^2)^2) , [0 Ls]);

    disp('j0...')
    j0 = chebfun(@(kx) besselj(0.0, kx), [0 Lo*Ls]);

    disp('j1...')
    j1 = chebfun(@(kx) besselj(1.0, kx), [0 max(Lo*Lx, Lo*Ls)]);

    disp('cos...')
    Cos = chebfun(@(lk) cos(lk), [0 Lo*Ll]);

    disp('sin...')
    Sin = chebfun(@(lk) sin(lk), [0 Lo*Ll]);

    xi = chebfun('x', [0 Lx]);

    omega = chebfun('x', [0 Lo]);

    disp('2...');

    inner_e_sl_2 = chebfun2(@(lambda, sigma) sum( j1(omega*sigma)/sigma * j1(omega*xi) * Sin(omega*lambda) ), [0 Ll 0.001 Ls], 'vectorize');

    figure(3)
    plot(inner_e_sl_2)

    disp('4...');

    u_sl = chebfun2( @(sigma, lambda) sum( xi^2*Phi(sigma) * inner_e_sl_2(lambda, sigma)), [0.001 Ls 0 Ll], 'vectorize');

    figure(4)
    plot(u_sl)

    disp('1...');

    inner_e_sl_1 = chebfun2(@(lambda, sigma) sum( j0(omega*sigma) * j1(omega*xi) * Cos(omega*lambda) ), [0 Ll 0.001 Ls], 'vectorize');


    disp('3...');

    eta_sl_1 = chebfun2( @(sigma, lambda) sum( xi^2*Phi(sigma) * inner_e_sl_1(lambda, sigma)), [0.001 Ls 0 Ll], 'vectorize');

    eta_sl = chebfun2(@(sigma,lambda) -1/4*eta_sl_1(sigma, lambda) -1/2*u_sl(sigma, lambda).^2, [0.001 Ls 0 Ll]);

    %eta_sl(1,1)

    figure(5)
    plot(eta_sl)

    %save('catalinaadhd');

end

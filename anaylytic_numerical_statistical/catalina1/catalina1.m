%warning: look at the parameters to make sure all are good

function [eta, u] = catalina1(save)

    global H1 H2 c1 c2 x1 x2 eta_0 eta_prime u_0 u_prime td t0 Tf x0 Xf

    %integration Parameters
    K = 18;
    La = 20;
    Ls = 20;

    %linear transform of x1, x2
    sigma1 = sqrt(x1*16);
    sigma2 = sqrt(x2*16);

    %eta(sigma, 0)
    eta_sigma_0 = @(sigma) H1*exp( -(c1*(sigma^2 - sigma1^2).^2)/256 ) - H2*exp( -(c2*(sigma^2 - sigma2^2).^2)/256 );

    eta_sigma_prime_0 = @(sigma) diff(eta_sigma_0(sigma));

    Phi = @(sigma) 4*eta_sigma_prime_0(sigma)/sigma;

    disp('j0...')
    j0 = chebfun(@(kx) besselj(0.0, kx), [0 K]);

    disp('j1...')
    j1 = chebfun(@(kx) besselj(1.0, kx), [0 K*10]);

    disp('cos...')
    Cos = chebfun(@(lk) cos(lk), [0 K*10], 'vectorize');

    disp('sin...')
    Sin = chebfun(@(lk) sin(lk), [0 K*10], 'vectorize');

    xi = chebfun('x', [0 K]);

    omega = chebfun('x', [0 K]);

    disp('2...');

    inner_e_sl_2 = chebfun3(@(lambda, sigma, xii) sum( j1(omega*sigma)/sigma * j1(omega*xii) * Sin(omega*lambda) ), [0 .1 0.01 .1 0 .1], 'vectorize');


    disp('1...');

    inner_e_sl_1 = chebfun3(@(lambda, sigma, xii) sum( j0(omega*sigma) * j1(omega*xii) * Cos(omega*lambda) ), [0 .1 0.01 .1 0 1.], 'vectorize');


    disp('3...');

    eta_sl_1 = chebfun2( @(sigma, lambda) sum( xi^2*Phi(sigma) * inner_e_sl_1(lambda, sigma, xi)), [0.01 .1 0 .1], 'vectorize');

    disp('4...');

    eta_sl_2 = chebfun2( @(sigma, lambda) sum( xi^2*Phi(sigma) * inner_e_sl_2(lambda, sigma, xi)), [0.01 .1 0 .1], 'vectorize');

    eta_sl = @(sigma,lambda) -1/4*eta_sl_1(sigma, lambda) -1/2*eta_sl_2(sigma, lambda)^2;

    eta_sl(1,1)

    if save
        disp('saving ...')
        save('eta_sl', 'eta_sl')
        disp('done')
    end
end

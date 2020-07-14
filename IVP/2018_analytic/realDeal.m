
%compute Nicolcsky 2018 anaylytic solution

function [psi, phi] = realDeal(data_proj, zero_inital_u, save_name)

    global eta_0 eta_prime u_0 u_prime t0 Tf x0 Xf upS lowS upL lowL K

    m = Inf;
    beta = 1;

    xx   = chebfun('x', [x0 25]);

    disp('j0...')
    j0 = chebfun(@(kx) besselj(0.0, kx), [0.0 max([2.0*K*max(sqrt(xx + eta_0(xx))) 2.0*K*sqrt(upS)])]);

    disp('j1...')
    j1 = chebfun(@(kx) besselj(1.0, kx), [0.0 max([2.0*K*max(sqrt(xx + eta_0(xx))) 2.0*K*sqrt(upS)])]);

    disp('cos...')
    Cos = chebfun(@(lk) cos(lk), [0 upL*K], 'vectorize');

    disp('sin...')
    Sin = chebfun(@(lk) sin(lk), [0 upL*K], 'vectorize');

    if data_proj

        disp('data projection to lambda = 0');

        disp('p...')
        x = chebfun('x', [0 upS]);

        s = @(x) x + eta_0(x);
        A = @(x) [0 1; beta^2*s(x) 0];
        B = [0 0; 1 0];

        D = @(x) (1+eta_prime(x))*eye(2) + u_prime(x)*A(x);

        phi0 = @(x) [u_0(x) ; eta_0(x)+(u_0(x).^2)/2];
        phi0_prime = @(x) [u_prime(x); eta_prime(x)+u_0(x)*u_prime(x)];

        proj = @(x) phi0(x) + u_0(x).*(u_prime(x).*A(x)*inv(D(x))*B*phi0(x) - B*phi0(x) - A(x)*inv(D(x))*phi0_prime(x));

        detD = @(x) (1+eta_0(x)).^2 - beta^2*s(x)*u_prime(x).^2;

        in1_phi = @(x) u_prime(x)*(1+eta_prime(x)*u_0(x));
        in2_phi = @(x) 0;
        in3_phi = @(x) -beta^2*s(x)*u_prime(x).^2 + (1+eta_prime(x))*(eta_prime(x)+u_0(x)*u_prime(x));
        proj_phi = chebfun(@(x) u_0(x) + u_0(x)*((in1_phi(x) - in3_phi(x))/detD(x) - in2_phi(x) ), [0 upS]);

        in1_psi = @(x) u_prime(x)*(-beta^.2*s(x)*u_0(x)^.2);
        in2_psi = @(x) u_0(x);
        in3_psi = @(x) beta^2*s(x)*(1+eta_prime(x))*u_prime(x) - u_prime(x)*beta^2*s(x)*(eta_prime(x)+u_0(x)*u_prime(x));
        proj_psi = chebfun(@(x) eta_0(x) + (u_0(x).^2)/2 + u_0(x)*((in1_psi(x) - in3_psi(x))/detD(x) - in2_psi(x) ), [0 upS]);

        %tests that show that the non-matrix transform produces the same
        %results as the matrix one
%         figure(1);
%         plot(proj_psi);
%
%         figure(2);
%         plot(proj_phi);
%
%         figure(3);
%         test1 = chebfun(@(k) psi_0(k), [0 Ls]);
%         plot(test1);
%
%         figure(4);
%         test2 = chebfun(@(k) phi_0(k), [0 Ls]);
%         plot(test2);


        disp('b...')

        b  = chebfun(@(k) -2*beta*k*sum( proj_phi(x)*(x+eta_0(x))^(1/2).*j1(2.0*k*sqrt(x + eta_0(x))).*(1 + eta_prime(x))), [0 K]);

        disp('a...');

        a = chebfun(@(k) 2*k*sum(proj_psi(x).*j0(2.0*k*sqrt(x + eta_0(x))).*(1 + eta_prime(x))), [0 K]);

    elseif zero_inital_u

        x = chebfun('x', [0 Xf]);

        disp('zero initial velocity so setting a....')
        a = chebfun(@(k) 2*k*sum(eta_0(x).*j0(2.0*k*sqrt(x + eta_0(x))).*(1 + eta_prime(x))), [0 K]);

        disp('b = 0....')
        b = chebfun(@(k) 0, [0 K]);

    else %use Phi 0

        disp('non zero_velocity and no data projection')
        disp('phi_0 and psi_0 ...')
        reg_phi_0 = @(x) u_0(x);
        reg_psi_0 = @(x) eta_0(x)+(u_0(x).^2)/2;

        p = chebfun('x', [0 K]);

        disp('a...')
        a  = chebfun(@(k) 2*k*sum( reg_psi_0(p)*j0(2*k*sqrt(p)) ), [0 K]);

        disp('b...')
        b  = chebfun(@(k) -2*beta*k*sum( reg_phi_0(p)*p^(1/2)*j1(2*k*sqrt(p)) ), [0 K]);

    end

    %phi and psi integration variable
    k = chebfun('x', [0 K]);

    disp('psi...')
    psi = chebfun2(@(s,la) sum( ( a(k)*Cos(la*k) + b(k)*Sin(la*k) ) * j0(2.0*k*sqrt(s)) ), [0 upS 0 upL], 'vectorize');

    disp('phi...')
    phi = chebfun2(@(s,la) s^(-1/2)*sum( ( a(k)*Sin(la*k) - b(k)*Cos(la*k) ) * j1(2.0*k*sqrt(s)) ), [0.00001 upS 0.00 upL], 'vectorize');

    disp('saving ...')
    save(save_name)
    disp('done')

    %phi_0 via matrix algerbra
    function phi = phi_0(x)
        phi = [1 0]*proj(x);
    end

    %psi_0 via matrix algerbra
    function psi = psi_0(x)
        psi = [0 1]*proj(x);
    end

end

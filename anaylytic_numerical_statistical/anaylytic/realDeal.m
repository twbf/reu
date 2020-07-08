
%compute Nicolcsky 2018 anaylytic solution

function [psi, phi] = realDeal(data_proj, zero_inital_u)

    %addpath('export_fig/');

    global H1 H2 c1 c2 x1 x2 eta_0 eta_prime u_0 u_prime t0 Tf x0 Xf

    %%% Physical parameter:
    K  = 1.0;		% upper bound in k domain [0, K]
    L  = 20.0;		% upper bound in x domain [0, L]
    Ls = Xf;		% upper bound for s parameter
    La = Tf;		% upper bound for lambda parameter

    %%% SWE Parameters

    m = Inf;
    beta = 1;

    %Plotting gamma
    figure(1);
    t = linspace(t0,Tf,100);
    fir = t+eta_0(t);
    sec = -u_0(t);
    %scatter(fir,sec)
    %export_fig('gamma.png', '-m2', '-a4', '-painters');

    xx   = chebfun('x', [0 L]);

    disp('j0...')
    j0 = chebfun(@(kx) besselj(0.0, kx), [0.0 max([2.0*K*max(sqrt(xx + eta_0(xx))) 2.0*K*sqrt(Ls)])]);

    disp('j1...')
    j1 = chebfun(@(kx) besselj(1.0, kx), [0.0 max([2.0*K*max(sqrt(xx + eta_0(xx))) 2.0*K*sqrt(Ls)])]);

    disp('cos...')
    Cos = chebfun(@(lk) cos(lk), [0 La*K], 'vectorize');

    disp('sin...')
    Sin = chebfun(@(lk) sin(lk), [0 La*K], 'vectorize');

    if data_proj

        disp('data projection to lambda = 0');

        disp('p...')
        x = chebfun('x', [0 Ls]);

        s = @(x) x + eta_0(x);
        A = @(x) [0 1; beta^2*s(x) 0];
        B = [0 0; 1 0];

        D = @(x) (1+eta_prime(x))*eye(2) + u_prime(x)*A(x);

        phi0 = @(x) [u_0(x) ; eta_0(x)+(u_0(x).^2)/2];
        phi0_prime = @(x) [u_prime(x); eta_prime(x)+u_0(x)*u_prime(x)];

        %proj = @(x) phi0(x) + u_0(x)*(u_prime(x)*A(x)*inv(D(x))*B*phi0(x) - B*phi0(x) - A(x)*inv(D(x))*phi0_prime(x));
        proj = @(x) phi0(x) + u_0(x).*(u_prime(x).*A(x)*inv(D(x))*B*phi0(x) - B*phi0(x) - A(x)*inv(D(x))*phi0_prime(x));

        detD = @(x) (1+eta_0(x)).^2 - beta^2*s(x)*u_prime(x).^2;

        in1_phi = @(x) u_prime(x)*(1+eta_prime(x)*u_0(x));
        in2_phi = @(x) 0;
        in3_phi = @(x) -beta^2*s(x)*u_prime(x).^2 + (1+eta_prime(x))*(eta_prime(x)+u_0(x)*u_prime(x));
        proj_phi = chebfun(@(x) u_0(x) + u_0(x)*((in1_phi(x) - in3_phi(x))/detD(x) - in2_phi(x) ), [0 Ls]);

        in1_psi = @(x) u_prime(x)*(-beta^.2*s(x)*u_0(x)^.2);
        in2_psi = @(x) u_0(x);
        in3_psi = @(x) beta^2*s(x)*(1+eta_prime(x))*u_prime(x) - u_prime(x)*beta^2*s(x)*(eta_prime(x)+u_0(x)*u_prime(x));
        proj_psi = chebfun(@(x) eta_0(x) + (u_0(x).^2)/2 + u_0(x)*((in1_psi(x) - in3_psi(x))/detD(x) - in2_psi(x) ), [0 Ls]);

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

        x = chebfun('x', [0 K]);

        disp('zero initial velocity so setting a....')
        a = chebfun(@(k) 2*k*sum(eta_0(x).*j0(2.0*k*sqrt(x + eta_0(x))).*(1 + eta_prime(x))), [0 K]);
        figure(2);
        plot(a);
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

    %plotting a and b

    %   figure(2);
    %   plot(a, '-', 'LineWidth', 2.0), grid off
    %   xlabel('$k$', 'interpreter', 'LaTeX', 'fontsize', 12);
    %   ylabel('$a(k)$', 'interpreter', 'LaTeX', 'fontsize', 12);
    %   title('Coefficient $a(k)$', 'interpreter', 'LaTeX', 'fontsize', 13);
    %   set(gcf, 'color', 'w');
    %   export_fig('a(k).png', '-m2', '-a4', '-painters');
    %
    %   figure(3);
    %   plot(b, '-', 'LineWidth', 2.0), grid off
    %   xlabel('$k$', 'interpreter', 'LaTeX', 'fontsize', 12);
    %   ylabel('$b(k)$', 'interpreter', 'LaTeX', 'fontsize', 12);
    %   title('Coefficient $b(k)$', 'interpreter', 'LaTeX', 'fontsize', 13);
    %   set(gcf, 'color', 'w');
    %   export_fig('b(k).png', '-m2', '-a4', '-painters');

    %phi and psi integration variable
    k   = chebfun('x', [0 K]);

    disp('psi...')
    psi   = chebfun2(@(s,la) sum( ( a(k)*Cos(la*k) + b(k)*Sin(la*k) ) * j0(2.0*k*sqrt(s)) ), [0 Ls 0 La], 'vectorize');

       figure(4);
       plot(psi);
    %   xlabel('$s$', 'interpreter', 'LaTeX', 'fontsize', 12);
    %   ylabel('$\lambda$', 'interpreter', 'LaTeX', 'fontsize', 12);
    %   view([0 90]); colorbar;
    %   title('Psi Two-parameters integral $f(s,\lambda)$', 'interpreter', 'LaTeX', 'fontsize', 12);
    %   set(gcf, 'color', 'w');
    %   export_fig('psi.png', '-m2', '-a4', '-painters');

    disp('phi...')
    phi = chebfun2(@(s,la) s^(-1/2)*sum( ( a(k)*Sin(la*k) - b(k)*Cos(la*k) ) * j1(2.0*k*sqrt(s)) ), [0.00001 Ls 0.00 La], 'vectorize');

       figure(5);
       plot(phi);
    %   xlabel('$s$', 'interpreter', 'LaTeX', 'fontsize', 12);
    %   ylabel('$\lambda$', 'interpreter', 'LaTeX', 'fontsize', 12);
    %   view([0 90]); colorbar;
    %   title('Phi Two-parameters integral $f(s,\lambda)$', 'interpreter', 'LaTeX', 'fontsize', 12);
    %   set(gcf, 'color', 'w');
    %   export_fig('phi.png', '-m2', '-a4', '-painters');

    disp('saving ...')
    save('psi_phi_projection_test')
    disp('done')

    function phi = phi_0(x)

        % input array x is a list of values for the projection
        % note: memory should be reaalocated

        phi = [1 0]*proj(x);

    end

    function psi = psi_0(x)

        % input array x is a list of values for the projection
        % note: memory should be reaalocated

        %for i=1:size(x,2)
            %tmp = proj(x(i));
           % psi(i) = tmp(2);
        %end

        psi = [0 1]*proj(x);

    end

end

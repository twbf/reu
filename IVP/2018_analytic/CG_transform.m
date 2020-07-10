% Goal: For a given phi(s, lambda) and psi(s, lambda) compute eta, u, x, t

function [ana_eta, ana_u] = CG_transform(loadd, size, data_proj, zero_inital_u)

    global H1 H2 c1 c2 x1 x2 eta_0 eta_prime u_0 u_prime td t0 Tf x0 Xf

    disp(" ")
    disp("Nicolcsky Anaylytic:");

    if loadd
        load psi_phi_projection_test2.mat phi psi
    else
        [psi, phi] = realDeal(data_proj, zero_inital_u);
    end


    %Parameters -----
    lambda_list = linspace(t0, Tf, size); %lambda space
    s = linspace(x0, Xf, size); %s space   note: at s = 0 regularization is needed - more on this later

    %deminsional variables
    g = 9.81; %gravity
    alpha = td; %slope
    l = 1; %arbitrary scaling parameter

    x = zeros(size, size);
    t = zeros(size, size);
    h = zeros(size, size);
    u = zeros(size, size);

    for i=1:size
        lambda = lambda_list(i); %for a single lambda

        %CG transform from (s, lambda) to (x,t)
        u(i,:) = real( phi(s,lambda) );
        h(i,:) = real( psi(s,lambda) )-u(i,:).^2/2;
        x(i,:) = s - h(i,:);
        t(i,:) = u(i,:) + lambda;

        %deminsionalizing
        u(i,:) = u(i,:)*sqrt(g*alpha*l);
        h(i,:) = h(i,:)*l*alpha;
        x(i,:) = x(i,:)*l;
        t(i,:) = t(i,:)/sqrt(g*alpha/l);
    end

    % scatter() and scatteredInterpolant() takes lists not matricies so x and t have to be reshaped
    xx = reshape(x, [size*size,1]);
    tt = reshape(t, [size*size,1]);
    hh = reshape(h, [size*size,1]);
    uu = reshape(u, [size*size,1]);

    ana_eta = scatteredInterpolant(xx,tt,hh);

    ana_u = scatteredInterpolant(xx,tt,uu);

    figure(1);
    scatter(xx,tt);
    %xlabel('$x$', 'interpreter', 'LaTeX', 'fontsize', 15);
    %ylabel('$t$', 'interpreter', 'LaTeX', 'fontsize', 15);
    %title('Equally Spaced Grid in $(s, \lambda)$ Transformed to $(x, t)$ ', 'interpreter', 'LaTeX', 'fontsize', 20);
    %export_fig('CG_coordinate_transform', '-m2', '-a4', '-painters');

    % figure(2);
    % mesh(x,t,h);
    % xlabel('$x$', 'interpreter', 'LaTeX', 'fontsize', 15);
    % ylabel('$t$', 'interpreter', 'LaTeX', 'fontsize', 15);
    % zlabel('$\eta$', 'interpreter', 'LaTeX', 'fontsize', 15);
    % title(' CG transformed $\eta$', 'interpreter', 'LaTeX', 'fontsize', 20);
    % export_fig('CG_eta', '-m2', '-a4', '-painters');
    %
    % figure(3);
    % mesh(x,t,u);
    % xlabel('$x$', 'interpreter', 'LaTeX', 'fontsize', 15);
    % ylabel('$t$', 'interpreter', 'LaTeX', 'fontsize', 15);
    % zlabel('$u$', 'interpreter', 'LaTeX', 'fontsize', 15);
    % title('CG transformed $u$', 'interpreter', 'LaTeX', 'fontsize', 20);
    % export_fig('CG_u', '-m2', '-a4', '-painters');
    %
    %
    % %interpolating mesh in order to get a smooth polynomial curve
    % analytic_eta = scatteredInterpolant(xx,tt,hh);
    % save('analytic_eta_s1', 'analytic_eta');
    %
    %
    % analytic_u = scatteredInterpolant(xx,tt,u);
    % save('analytic_u', 'analytic_u');

end

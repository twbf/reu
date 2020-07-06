function [eta, u] = catalina_transform(save, size)

    global H1 H2 c1 c2 x1 x2 eta_0 eta_prime u_0 u_prime td t0 Tf x0 Xf

    disp(" ")
    disp("Catalina1 Anaylytic:");

    %this returns eta and u which are functions of sigma and lambda
    [eta_sl, u_sl] = catalina1(save);

    %Parameters -----
    lambda_list = linspace(t0, Tf, size); %lambda space
    s = linspace(x0, Xf, size); %s space   note: at s = 0 regularization is needed - more on this later

    %deminsional variables
    g = 9.81; %gravity
    beta = tan(td); %slope
    l = 1; %arbitrary scaling parameter

    x = zeros(size, size);
    t = zeros(size, size);
    h = zeros(size, size);
    u = zeros(size, size);

    for i=1:size
        lambda = lambda_list(i); %for a single lambda

        %CG transform from (s, lambda) to (x,t)
        u(i,:) = u_sl(s,lambda);
        h(i,:) = eta_sl(s, lambda);
        x(i,:) = s - h(i,:);
        t(i,:) = u(i,:) + lambda;

        %deminsionalizing
        u(i,:) = u(i,:)*sqrt(g*beta*l);
        h(i,:) = h(i,:)*l*beta;
        x(i,:) = x(i,:)*l;
        t(i,:) = t(i,:)*sqrt(l/(beta*g));
    end

    % scatter() and scatteredInterpolant() takes lists not matricies so x and t have to be reshaped
    xx = reshape(x, [size*size,1]);
    tt = reshape(t, [size*size,1]);
    hh = reshape(h, [size*size,1]);
    uu = reshape(u, [size*size,1]);

    eta = scatteredInterpolant(xx,tt,hh);

    u = scatteredInterpolant(xx,tt,uu);

end

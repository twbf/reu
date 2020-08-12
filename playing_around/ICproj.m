% A program that computes U_n(x) through data projection

close all; clear all; clc;
format long


g  = 9.8;	 % gravity acceleration
d = 1;

for i=1:2
    % Projection parameters
    n = 4;
    ChebDiv = 5;
    chebN = i;

    % Spatial grid parameters
    a  = 0.0;	% the left boundary (incident wave)
    b  = 10.0;	% the right boundary (wall)
    Ndx = 2000; % spatial step
    x = linspace(a,b,Ndx);

    % Initial condition parameters (eta):
    H1 = 0.006;
    H2 = 0.018;
    c1 = 0.4444;
    c2 = 4.0;
    x1 = 4.1209;
    x2 = 1.6384;

    % Inital eta function
    eta0Dim = zeros(1,Ndx);
    eta0Dim(1:Ndx) = H1*exp(-c1*(x - x1).^2) - H2*exp(-c2*(x - x2).^2);

    % Bathymetry
    td = 10.0/10.0; % bottom slope
    h = td*x;

    % initial u function
    u0Dim = zeros(1,Ndx);
    u0Dim(1:Ndx) = sqrt(g*(eta0Dim+h));

    % removing dimensions
    eta0 = eta0Dim/(d);            % dimensionless eta
    u0 = u0Dim/(sqrt(g*d));        % dimensionless u

    % computing sigma(x) at t=0
    s0 = x+eta0;

    % fitting chebfun to sigma(x) and creating inverse function
    s0Cheb = polyfit(x(1,:),s0(1,:),chebN,domain(a,b));
    gs0Cheb = polyfit(s0(1,:),x(1,:),chebN,domain(a,b));    % gamma(sigma)
    gs0 = gs0Cheb(x(1,:));

    % fitting chebfun to u0 and eta0
    u0Cheb = polyfit(x(1,:),u0(1,:),chebN,domain(a,b));
    eta0Cheb = polyfit(x(1,:),eta0(1,:),chebN,domain(a,b));
    % evaluating u0(gamma(sigma))
    u0gs = u0Cheb(gs0);
    eta0gs = eta0Cheb(gs0);

    % fitting chebfun for u(gamma(sigma)) and eta(gamma(sigma))
    u0gsCheb = polyfit(s0(1,:),u0gs(1,:),chebN,domain(a,b));
    eta0gsCheb = polyfit(s0(1,:),eta0gs(1,:),chebN,domain(a,b));

    % components of phi vector
    Phi0TCheb = u0gsCheb;
    Phi0BCheb = eta0gsCheb+0.5*u0gsCheb.^2;
    Phi0T = Phi0TCheb(x);
    Phi0B = Phi0BCheb(x);

    % derivative of Phi0T
    dPhi0TdsCheb = diff(Phi0TCheb);
    dPhi0BdsCheb = diff(Phi0BCheb);
    dPhi0Tds = dPhi0TdsCheb(x);
    dPhi0Bds = dPhi0BdsCheb(x);

    % creating and inverting matrix D
    one = repelem(1,Ndx);
    invD = zeros(4,Ndx);
    for i=1:Ndx
        matrixD = [one(i),dPhi0Tds(i);s0(i)*dPhi0Tds(i),one(i)];
        invD = inv(matrixD);
        invLT(1,i) = invD(1,1);                   % left-top element
        invRT(1,i) = invD(1,2);                   % right-top element
        invLB(1,i) = invD(2,1);                   % left-bottom element
        invRB(1,i) = invD(2,2);                   % right-bottom element
        invD = [invLT;invRT;invLB;invRB];         % inverted matrix D (4xM matrix)
    end

    % empty arrays
    FkTsave = zeros(n,Ndx);         % store F_k top
    FkBsave = zeros(n,Ndx);         % store F_k bottom
    FkTds = zeros(1,Ndx);           % derivative top
    FkBds = zeros(1,Ndx);           % derivative bottom
    sumStoreT = zeros(n,Ndx);       % storing sum for each K value
    sumStoreB = zeros(n,Ndx);
    Phi_nT = zeros(1,Ndx);          % solution arrays
    Phi_nB = zeros(1,Ndx);

    % main loop
    for i=1:n
        term1 = ((Phi0T).^i)/factorial(i);
        FkTds = dPhi0Tds;
        FkBds = dPhi0Bds;
        kphi = i*Phi0T;
        nDinv_Fk = [FkTds;kphi+FkBds];
        FkT = invD(1,:).*nDinv_Fk(1,:)+invD(2,:).*nDinv_Fk(2,:);
        FkB = invD(3,:).*nDinv_Fk(1,:)+invD(4,:).*nDinv_Fk(2,:);
        FkTsave(i,:) = FkT;
        FkBsave(i,:) = FkB;
        dPhi0TdsCheb = polyfit(s0(1,:),FkT(1,:),chebN,domain(a,b));
        dPhi0BdsCheb = polyfit(s0(1,:),FkB(1,:),chebN,domain(a,b));
        dPhi0Tds = dPhi0TdsCheb(x);
        dPhi0Bds = dPhi0BdsCheb(x);
        sumStoreT(i,:) = term1.*FkT;
        sumStoreB(i,:) = term1.*FkB;
    end


    %Phi_nT = smooth(Phi0T+sum(sumStoreT));
    %Phi_nB = smooth(Phi0B+sum(sumStoreB));


    figure(1)
    plot(s0, Phi_nT,'.-','LineWidth',1)


    figure(2)
    plot(s0, Phi_nB,'.-','LineWidth',1)
    pause(1)
end


function Phi = order_n_dp(x, n)

    % Global variables:
    global eta_0 u_0
    global x_res t_res Xf g td

    eta0 = zeros(1,x_res);
    eta0(1:x_res) = eta_0(x);

    u0 = zeros(1,x_res);
    u0(1:x_res) = u_0(x);

    % computing sigma(x) at t=0
    s0 = x+eta0;

    % interpolating
    gs0Cheb = chebfun.interp1(s0(1,:),x(1,:), 'pchip');
    u0Cheb = chebfun.interp1(x(1,:),u0(1,:),'pchip');
    eta0Cheb = chebfun.interp1(x(1,:),eta0(1,:),'pchip');
    gs0 = gs0Cheb(x(1,:));

    % evaluating u0(gamma(sigma)) and eta0(gamma(sigma))
    u0gs = u0Cheb(gs0);
    eta0gs = eta0Cheb(gs0);

    eta0gsCheb = chebfun.interp1(s0(1,:),eta0gs(1,:),'pchip');
    u0gsCheb = chebfun.interp1(s0(1,:),u0gs(1,:),'pchip');

    Phi0TCheb = u0gsCheb;
    Phi0BCheb = eta0gsCheb+0.5*u0gsCheb.^2;
    Phi0T = Phi0TCheb(gs0);
    Phi0B = Phi0BCheb(gs0);

    % derivative of Phi0T
    dPhi0TdsCheb = diff(Phi0TCheb);
    dPhi0BdsCheb = diff(Phi0BCheb);
    dPhi0Tds = dPhi0TdsCheb(gs0);
    dPhi0Bds = dPhi0BdsCheb(gs0);


    % creating and inverting matrix D
    one = repelem(1,x_res);
    invD = zeros(4,x_res);
    for i=1:x_res
      matrixD = [one(i),dPhi0Tds(i);s0(i)*dPhi0Tds(i),one(i)];
      invD = inv(matrixD);
      invLT(1,i) = invD(1,1);                   % left-top element
      invRT(1,i) = invD(1,2);                   % right-top element
      invLB(1,i) = invD(2,1);                   % left-bottom element
      invRB(1,i) = invD(2,2);                   % right-bottom element
      invD = [invLT;invRT;invLB;invRB];         % inverted matrix D (4xM matrix)
    end

    % empty arrays
    FkTsave = zeros(n,x_res);         % store F_k top
    FkBsave = zeros(n,x_res);         % store F_k bottom
    FkTds = zeros(1,x_res);           % derivative top
    FkBds = zeros(1,x_res);           % derivative bottom
    sumStoreT = zeros(n,x_res);       % storing sum for each K value
    sumStoreB = zeros(n,x_res);
    Phi_nT = zeros(1,x_res);          % solution arrays
    Phi_nB = zeros(1,x_res);

    % initializing FkT
    FkT = Phi0T;

    % main loop
    for i=1:n
      term1 = ((Phi0T).^i)/factorial(i);
      FkTds = dPhi0Tds;
      FkBds = dPhi0Bds;
      kphi = i*FkT;
      nDinv_Fk = [FkTds;kphi+FkBds];
      FkT = invD(1,:).*nDinv_Fk(1,:)+invD(2,:).*nDinv_Fk(2,:);
      FkB = invD(3,:).*nDinv_Fk(1,:)+invD(4,:).*nDinv_Fk(2,:);
      FkTsave(i,:) = FkT;
      FkBsave(i,:) = FkB;
      dPhi0TdsCheb = chebfun.interp1(s0(1,:),FkT(1,:),'pchip');
      dPhi0BdsCheb = chebfun.interp1(s0(1,:),FkB(1,:),'pchip');
      dPhi0Tds = dPhi0TdsCheb(gs0);
      dPhi0Bds = dPhi0BdsCheb(gs0);
      sumStoreT(i,:) = term1.*FkT;
      sumStoreB(i,:) = term1.*FkB;
    end

    Phi_nT = Phi0T+sum(sumStoreT);
    Phi_nB = Phi0B+sum(sumStoreB);

    Phi = [Phi_nT; Phi_nB];

end

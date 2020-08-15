
function Psi = order_n_BC_proj(eta1, u1, n)

  global g t0 Tf t_res

  % parameters
  d=1;

  u1_dim = u1;
  eta1_dim = eta1;

  % removing dimensions
  t_dim = linspace(t0,Tf,t_res);
  t = t_dim/(sqrt(d/g));          % dimensionless time

  eta = eta1_dim/(d);             % dimensionless eta
  u = u1_dim/(sqrt(g*d));         % dimensionless u

  % computing lambda
  lambda = t-u;                % lambda array

  % polynomial fitting to find gamma(lambda)
  gammaCheb = chebfun.interp1(lambda,t,'pchip');
  gamma = gammaCheb(lambda);

  % computing eta at gamma
  etaGammaCheb = chebfun.interp1(gamma,eta,'pchip');
  etaGamma = etaGammaCheb(lambda);

  sigmaCheb = etaGammaCheb+1;

  Psi0T = gamma-lambda;                       % top component (phi)
  Psi0B = etaGamma+0.5*(gamma-lambda).^2;     % bottom component (psi)
  Psi0 = [Psi0T;Psi0B];

  Psi0TCheb = chebfun.interp1(lambda,Psi0T,'pchip');
  Psi0BCheb = chebfun.interp1(lambda,Psi0B,'pchip');
  dPsi0TdsCheb = diff(Psi0TCheb);
  dPsi0BdsCheb = diff(Psi0BCheb);
  dPsi0Tds = dPsi0TdsCheb(lambda);
  dPsi0Bds = dPsi0BdsCheb(lambda);

  dsdlCheb = diff(sigmaCheb);
  dsdl = dsdlCheb(lambda);

  sl = sigmaCheb(lambda);
  minusone = repelem(-1,t_res);

  % creating and inverting matrix D
  invD = zeros(4,t_res);
  for i=1:t_res
      matrixD = [dsdl(i),minusone(i);-sl(i),dsdl(i)];
      invD = inv(matrixD);
      invLT(1,i) = invD(1,1);                   % left-top element
      invRT(1,i) = invD(1,2);                   % right-top element
      invLB(1,i) = invD(2,1);                   % left-bottom element
      invRB(1,i) = invD(2,2);                   % right-bottom element
      invD = [invLT;invRT;invLB;invRB];         % inverted matrix D (4xM matrix

      % checking the deteriminant
      det0 = det(matrixD);
      if det0 == 0
          disp('detD=0')
      end
  end

  % empty arrays
  FkTsave = zeros(n,t_res);         % store F_k top
  FkBsave = zeros(n,t_res);         % store F_k bottom
  FkTds = zeros(1,t_res);           % derivative top
  FkBds = zeros(1,t_res);           % derivative bottom
  sumStoreT = zeros(n,t_res);       % storing sum for each K value
  sumStoreB = zeros(n,t_res);
  Psi_nT = zeros(1,t_res);          % solution arrays
  Phi_nB = zeros(1,t_res);

  % initializing FkT
  FkT = Psi0T;

  % main loop
  for i=1:n
      term1 = ((1-sl).^i)/factorial(i);
      FkTds = dPsi0Tds;
      FkBds = dPsi0Bds;
      kphi = i*FkT;
      nDinv_Fk = [FkTds;kphi+FkBds];
      FkT = invD(1,:).*nDinv_Fk(1,:)+invD(2,:).*nDinv_Fk(2,:);
      FkB = invD(3,:).*nDinv_Fk(1,:)+invD(4,:).*nDinv_Fk(2,:);
      FkTsave(i,:) = FkT;
      FkBsave(i,:) = FkB;
      dPsi0TdsCheb = chebfun.interp1(lambda(1,:),FkT(1,:),'pchip');
      dPsi0BdsCheb = chebfun.interp1(lambda(1,:),FkB(1,:),'pchip');
      dPhi0Tds = dPsi0TdsCheb(lambda);
      dPhi0Bds = dPsi0BdsCheb(lambda);
      sumStoreT(i,:) = term1.*FkT;
      sumStoreB(i,:) = term1.*FkB;
  end

  Psi_nT = Psi0T+sum(sumStoreT);
  Psi_nB = Psi0B+sum(sumStoreB);

  Psi = [Psi_nT; Psi_nB];

end

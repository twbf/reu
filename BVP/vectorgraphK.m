%--------- Graphing Boundary Data Vector---------%

close all
clear all
format Long

load('str_sine.mat','shelf')

%----Paramters------%
g = 9.81;
g2 = 0.5*g;
d=1;

t0 = 0;
tf = 10;

k = 10;
%---------Removing Dimensions------%
u1_dim = shelf.u1;                      % pulls u1 from runwave.m
u1 = u1_dim/(sqrt(g*d));                % dimensionless velocity array

Ndt = length(u1);                       % number of time steps
t_dim = linspace(t0, tf, Ndt);          % time array
t = t_dim/(sqrt(d/g));                  % dimensionless time array

eta1_dim = shelf.eta1;                  % pulls eta from runwave.m
eta1 = eta1_dim/(d);                    % making height dimensionless

n = length(u1);                         % number of array elements
%-------lambda-------%
lambda = t-u1;

%-------Polynomial Fitting to find Gamma-------%
coeff = polyfit(lambda,t,1);            % fitting a function to the inverse
c1 = coeff(1);
c2 = coeff(2);

gamma = c1*lambda+c2;  % gamma vector
gamma(1) = 0;
gamma(n) = gamma(n-1);

%-----Calculating Eta at Gamma------%
eta_of_gamma = interp1(t,eta1,gamma);   % calculating eta at gamma

%---------Sigma at boundary---------%
sigma = eta_of_gamma+1;                 % finding sigma

%---------computing psi vector------%
Psi0T = gamma-lambda;
Psi0B = eta_of_gamma+0.5*(gamma-lambda).^2;
Psi0 = [Psi0T; Psi0B];

%--------Kth derivative of sigma------------%

%{
Stores sigma as first row of matrix. Kth order derivatives are stored
in rows below.
%}

dsdl = zeros(k+1,n);
dsdl(1,:) = sigma;
for i=1:k+1
  dsdl(i+1,1) = (dsdl(i,2)-dsdl(i,1))/(lambda(2)-lambda(1));    % initial point
  for j=2:n-1
    dsdl(i+1,j) = (dsdl(i,j+1)-dsdl(i,j-1))/(lambda(j+1)-lambda(j-1));
  end
  dsdl(i+1,n) = (dsdl(i,n)-dsdl(i,n-1))/(lambda(n)-lambda(n-1)); % final point
end

%------Inverting Matrix D------%
minusone = repelem(-1, n);

%{
Takes Matrix D (2x2) and inverts the elements. The inverted elements
are then stored in a 4x200 matrix where row 1-4 correspond to elements
1-4 of matrix D.
%}
Dinverse = zeros(4,n);

for i=1:n
  Dmatrix = [dsdl(2,i),minusone(i);-sigma(i),dsdl(2,i)];
  Dmatrix_inv = inv(Dmatrix);

  LTinv(1,i) = Dmatrix_inv(1,1);         %LT
  RTinv(1,i) = Dmatrix_inv(1,2);         %RT
  LBinv(1,i) = Dmatrix_inv(2,1);         %LB
  RBinv(1,i) = Dmatrix_inv(2,2);         %RB
  Dinv = [LTinv; RTinv; LBinv; RBinv];
end
    % D inverse as 4x200 matrix
%--------Psi_n of lambda-------%

% Empty Arrays
FkTsave = zeros(k,n);         % store F_k top
FkBsave = zeros(k,n);         % store F_k bottom

FkTdl = zeros(2,n);       % derivative top
FkBdl = zeros(2,n);       % derivative bottom

FkTdl(1,:) = Psi0T;       % initializing derivative arrays
FkBdl(1,:) = Psi0B;

sumStoreT = zeros(k,n);   % storing sum for each K value
sumStoreB = zeros(k,n);

Psi_nT = zeros(1,n);      % solution arrays
Psi_nB = zeros(1,n);

for j=1:k

  % FIRST TERM IN SUMATION
  term1_fac = factorial(j);
  term1 = ((1-sigma).^j)/term1_fac;

  % COMPUTING DERIVATIVES
  for T=1:2
    FkTdl(T+1,1) = (FkTdl(T,2)-FkTdl(T,1))/(lambda(2)-lambda(1));    % initial point
    for Tdl=2:n-1
      FkTdl(T+1,Tdl) = (FkTdl(T,Tdl+1)-FkTdl(T,Tdl-1))/(lambda(Tdl+1)-lambda(Tdl-1));
    end
    FkTdl(T+1,n) = (FkTdl(T,n)-FkTdl(T,n-1))/(lambda(n)-lambda(n-1)); % final point
  end

  for B=1:2
    FkBdl(B+1,1) = (FkBdl(B,2)-FkBdl(B,1))/(lambda(2)-lambda(1));    % initial point
    for Bdl=2:n-1
      FkBdl(B+1,Bdl) = (FkBdl(B,Bdl+1)-FkBdl(B,Bdl-1))/(lambda(Bdl+1)-lambda(Bdl-1));
    end
    FkBdl(B+1,n) = (FkBdl(B,n)-FkBdl(B,n-1))/(lambda(n)-lambda(n-1)); % final point
  end

  % K TIME TOP ELEMENT OF PSI
  KxPhi = j*Psi0T;

  % COMPUTING PART OF F_K. STILL NEED TO MULTIPLY BY D INVERSE.
  nDinv_Fk = [FkTdl(2,:); KxPhi+FkBdl(2,:)];

  % MULTIPLYING BY D INVERSE
  FkT = Dinv(1,:).*nDinv_Fk(1,:)+Dinv(2,:).*nDinv_Fk(2,:);
  FkB = Dinv(3,:).*nDinv_Fk(1,:)+Dinv(4,:).*nDinv_Fk(2,:);

  FkTsave(j,:) = FkT;
  FkBsave(j,:) = FkB;

  % UPDATING VECTOR OPERATED UPON TO F_K-1
  FkTdl(1,:) = FkT;
  FkBdl(1,:) = FkB;

  % SUMATION FOR EACH K
  sumStoreT(j,:) = term1.*FkT;
  sumStoreB(j,:) = term1.*FkB;

end

% COMPUTING PSI_N
for i=1:n
  Psi_nT(i) = Psi0T(i)+sum(sumStoreT(i));
  Psi_nB(i) = Psi0B(i)+sum(sumStoreB(i));
end

Psi_n = [Psi_nT;Psi_nB];                  %Psi_n

sumexpT = sumStoreT;
sumexpB = sumStoreB;

for i=1:k
  figure(i)
  plot(lambda, Psi0B+sum(sumexpB(1:i,:)))
  axis([0 30 0.99 1.01])
  pause(0.5)
end

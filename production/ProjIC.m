close all; clear all; clc;
format long

% Global variables:
global cf2 d xc FS IN LW
global a amp td b g g2 dx h N

g  = 9.8;	 % gravity acceleration
d = 1;

n = 5;

% Spatial grid parameters
a  = 0.0;	% the left boundary (incident wave)
b  = 10.0;	% the right boundary (wall)
Ndx = 10000; % spatial step
x = linspace(a,b,Ndx);

% Initial condition parameters (eta):
H1 = 0.006;
H2 = 0.018;
c1 = 0.4444;
c2 = 4.0;
x1 = 4.1209;
x2 = 1.6384;

eta0Dim = zeros(1,Ndx);
eta0Dim(1:Ndx) = H1*exp(-c1*(x - x1).^2) - H2*exp(-c2*(x - x2).^2);

% Bathymetry
td = 10.0/10.0; % bottom slope
h = td*x;

% initial u function
u0Dim = zeros(1,Ndx);
u0Dim(1:Ndx) = 0.1*sin(5*x).*exp(-0.4*(x-5).^2);

% removing dimensions
eta0 = eta0Dim/(d);            % dimensionless eta
u0 = u0Dim/(sqrt(g*d));        % dimensionless u

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

figure(1)
plot(s0, Phi_nT,'.-','LineWidth',1)


figure(2)
plot(s0, Phi_nB,'.-','LineWidth',1)

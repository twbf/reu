%{
    Finite volume solver for NSWE wave propagation and run-up
    Copyright (C) 2020 Denys DUTYKH

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

%}

%%% -------------------------------------------------- %%%
%%% The main file for NSWE solver based on FV method   %%%
%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, CNRS -- LAMA, Univ of Savoie %%%
%%% E-mail: Denys.Dutykh@univ-savoie.fr                %%%
%%% Web:    http://www.denys-dutykh.com/               %%%
%%% Blog:   http://dutykh.github.io/                   %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%
%%% Distributed under GNU General Public License       %%%
%%% -------------------------------------------------- %%%

clc
clear all
close all
format longE

%%% Libraries we use:
addpath('sources/');
addpath('odetpbar/');
addpath('chebfun-master');

%%% Global variables:
global cf2 d xc FS IN LW
global a amp td b g g2 dx h N

%--------RESOLUTION PARAMETERS------------%

M  = 10000	% number of time instances where we project solution
N  = 3000   % number of grid points


%-------------DATA PROJ PARAMETERS-------------%

n = 0;   % order of data projection (k)
ChebDiv = 2; % number of divisions in Chebfun
chebN = 25; % number of cheb points


%--------------HODOGRAPH SOLVER PARAMETERS---------%

numSig = 15;              % number of sigma points
numLam = 10000;          % number of lambda points

sig = linspace(0,1,numSig);
lam = linspace(0,2,numLam); %lambda, leave at 0-10

dLam = lam(2)-lam(1);
dSig = sig(2)-sig(1);

courant = dLam/dSig

%------------------------------------------------------------

FS = 'FontSize';
IN = 'Interpreter';
LS = 'LineStyle';
LW = 'LineWidth';

%%% Physical parameters:
g  = 9.8;	% gravity acceleration
g2 = 0.5*g;	% g/2
cf2 = 0.0;	% friction coefficient

a  = -2.0;	% the left boundary (incident wave)
b  = 10.0;	% the right boundary (wall)
d = 1;
% Bottom slope:
td = 10.0/10.0;

%%% Initial condition parameters:

H1 = 0.006;
H2 = 0.018;
c1 = 0.4444;
c2 = 4.0;
x1 = 4.1209 +3;
x2 = 1.6384 +3;

%%% Numerical parameters:
x  = linspace(a, b, N+1)';			% cell interfaces
dx = x(2) - x(1);                  	% spatial grid step
xc = 0.5*(x(1:end-1) + x(2:end));  	% centers of cells

%%% Bathymetry function:
h  = td*xc;

%%% Choice of the initial condition:
w0 = zeros(2*N,1);
w0(1:N) = max(h, eps+0*h);   % zero initial condition without velocities
w0(1:N) = w0(1:N) + H1*exp(-c1*(xc - x1).^2) - H2*exp(-c2*(xc - x2).^2);

% time stepping:
t0 = 0.0;
Tf = 2.0; % final simulation time (the same as experiment)

%%% Plot the initial condition to check:
% amp = 1.5*H2;
% figure(5);
% set(gcf, 'pos', [1 621 903 353]);
% Plot(t0, w0);

%%% We run the simulation:
options = odeset('AbsTol', 1e-4, 'RelTol', 1e-4, 'OutputFcn', @odetpbar, 'MaxStep', 1.0);
fprintf('Simulation...\n');
sol = ode23(@RHS, [t0 Tf], w0, options);
fprintf(' Done\n');


%%% Post-processing of the solution:
tlist = linspace(t0, Tf, M);
solpr = deval(sol, tlist);


% %%% Make a little animation of the numerical solution:
Rup = zeros(M,1);

% figure(22)
% for t=1:M % loop in time
%     ind = find(solpr(1:N,t) > 1e-2, 1, 'first');
%     Rup(t) = -h(ind);
%
% 	Plot(tlist(t), solpr(:,t));
%     %disp(solpr(2501,t));                    % eta indexed at 501
%                                              % u indexed at 2501
%
% end % for t


xeq=x==1;
xeq1_eta = find(xeq) ;  %finds index
xeq1_u = find(xeq)+N;

u1_dim = solpr(xeq1_u,:);
eta1_dim = solpr(xeq1_eta,:);

% phi = solpr(1:M,1:M);
% psi = solpr(M+1:end,1:end);

xT = x.';
xSec = xT(1:end-1);

% figure(98)
% mesh(tlist,xSec,phi)
% title('U')
%
% figure(99)
% mesh(tlist,xSec,psi)
% title('Eta')
% xlabel('t')
% ylabel('psi')

Rup = Rup - Rup(1);

% figure (23)
%  for t=1:M % loop in time
%      ind1 = find(solpr(501,t) > 1e-2, 1, 'first');
%      eta(t) = -h(ind1);
%
%  	Plot(tlist(t), solpr(501,t));
%  end % for t



%%% Extraction of run-up data:
% Rup = smooth(Rup, 7);
% figure (24)
% plot(tlist, Rup, '-', LW, 2.0); grid off;
% axis tight;
% xlim([t0 Tf]); ylim([1.2*min(Rup) 1.1*max(Rup)]);
% xlabel('$t$, s', IN, 'LaTeX', FS, 14);
% ylabel('$R(t)$, m', IN, 'LaTeX', FS, 14);
% ll = legend(' Numerical data', 'location', 'SouthEast');
% set(ll, 'box', 'off');
% title('Wave run-up');
% set(gcf, 'Color', 'w');


% figure (25)
% plot(tlist, eta, '-', LW, 2.0); grid off;
% axis tight;
% xlim([t0 Tf]); ylim([1.2*min(Rup) 1.1*max(Rup)]);
% xlabel('$t$, s', IN, 'LaTeX', FS, 14);
% ylabel('$R(t)$, m', IN, 'LaTeX', FS, 14);
% ll = legend(' Numerical data', 'location', 'SouthEast');
% set(ll, 'box', 'off');
% title('BC');
% set(gcf, 'Color', 'w');

%--------------PSI_N_CHEBFUN START-------------------%

t_dim = linspace(t0,Tf,M);
t = t_dim/(sqrt(d/g));          % dimensionless time

eta = eta1_dim/(d);             % dimensionless eta

u = u1_dim/(sqrt(g*d));    %dimensionless u

figure(1)
plot(eta)

figure(2)
plot(u)



% computing lambda
lambda = t-u;                % lambda array

% polynomial fitting to find gamma(lambda)
gammaCoeff = polyfit(lambda, t, 1);
gamma = polyval(gammaCoeff, t);

% computing eta at gamma
etaGamma = interp1(t,eta,gamma);
etaGamma(end) = etaGamma(end-1);

% sigma at boundary
sigma = etaGamma+1;

% computing psi_0 vector
Psi0T = gamma-lambda;                     % top component (phi)
Psi0B = etaGamma+0.5*(gamma-lambda).^2;   % bottom component (psi)
Psi0 = [Psi0T;Psi0B];
% reshaping Psi0T,Psi0B,lambda and sigma for fitting multiple chebfuns

column = ChebDiv;                            % number of chebfuns to be fit
row = M/ChebDiv;                             % number of points being fit

while row ~= floor(row)
  column = column+1;                 % ensures number of columns is an integer
  row = M/column;
end

Psi0T = reshape(Psi0T, [row, column]);
Psi0B = reshape(Psi0B, [row, column]);
lambda = reshape(lambda, [row, column]);
sigma = reshape(sigma, [row, column]);

% Repeating elements
% Elements at tops of columns moved to bottom and bottoms moved to top to create
% overlap in domains and reduce error in further derivatives
domOL = 100;                          % domain overlap
rowOL = row+2.*domOL;               % row overlap

repT = 1:domOL;                     % identifies top elements to replace
repB = (rowOL-domOL+1):rowOL;       % identifies bottom elements to replace
idT = rowOL-2.*domOL+1:rowOL-domOL; % identifies top elements of original matrix
idB = domOL+1:2.*domOL;             % identifies bottom elements of original matrix

zB = zeros(domOL, column);          % added to bottom of matrix
zT = zeros(domOL, column);          % added to top of matrix

lambdaOL = [zT;lambda;zB];                  % modified matricies
sigmaOL = [zT;sigma;zB];
Psi0BOL = [zT;Psi0B;zB];
Psi0TOL = [zT;Psi0T;zB];

lambdaOL(repB,column) = lambdaOL(rowOL-domOL, column);   % assigning values to
sigmaOL(repB,column) = sigmaOL(rowOL-domOL, column);     % end of last column
Psi0BOL(repB,column) = Psi0BOL(rowOL-domOL, column);
Psi0TOL(repB,column) = Psi0TOL(rowOL-domOL, column);

lambdaOL(repT,1) = lambdaOL(domOL+1,1);             % assigning values to
sigmaOL(repT,1) = sigmaOL(domOL+1,1);               % begining of first column
Psi0BOL(repT,1) = Psi0BOL(domOL+1,1);
Psi0TOL(repT,1) = Psi0TOL(domOL+1,1);

for i=1:column-1
  lambdaOL(repB,i) = lambdaOL(idB,i+1);       % lambdaOL
  lambdaOL(repT,i+1) = lambdaOL(idT,i);
  sigmaOL(repB,i) = sigmaOL(idB,i+1);         % sigmaOL
  sigmaOL(repT,i+1) = sigmaOL(idT,i);
  Psi0TOL(repB,i) = Psi0TOL(idB,i+1);         % Psi0TOL
  Psi0TOL(repT,i+1) = Psi0TOL(idT,i);
  Psi0BOL(repB,i) = Psi0BOL(idB,i+1);         % Psi0BOL
  Psi0BOL(repT,i+1) = Psi0BOL(idT,i);
end

% fitting chebfuns to Psi0T, Psi0B, and sigma (as functions of lambda)
for i=1:column
  l0_OL = lambdaOL(1,i);              % begining of overlaping domain
  lf_OL = lambdaOL(rowOL,i);          % end of overlaping domain
  l0 = lambda(1,i);                   % begining of domain
  lf = lambda(row,i);                 % end of domain

  % overlaping chebfuns
  chebPsi0BOL = polyfit(lambdaOL(:,i),Psi0BOL(:,i),chebN,domain(l0_OL,lf_OL));
  chebPsi0TOL = polyfit(lambdaOL(:,i),Psi0TOL(:,i),chebN,domain(l0_OL,lf_OL));
  chebSigmaOL = polyfit(lambdaOL(:,i),sigmaOL(:,i),chebN,domain(l0_OL,lf_OL));

  % saving overlaping chebfuns
  save(sprintf('chebPsi0BOL_%d.mat',i), 'chebPsi0BOL');
  save(sprintf('chebPsi0TOL_%d.mat',i), 'chebPsi0TOL');
  save(sprintf('chebSigmaOL_%d.mat',i), 'chebSigmaOL');

  % restricting functions so domains do not overlap
  chebPsi0B = restrict(chebPsi0BOL, [l0, lf]);
  chebPsi0T = restrict(chebPsi0TOL, [l0, lf]);
  chebSigma = restrict(chebSigmaOL, [l0, lf]);

  % saving chebfuns
  save(sprintf('chebPsi0B_%d.mat',i), 'chebPsi0B');
  save(sprintf('chebPsi0T_%d.mat',i), 'chebPsi0T');
  save(sprintf('chebSigma_%d.mat',i), 'chebSigma');
end

% creating elements of matrix D

% derivative of sigma with respect to lambda (chebfun)
for i=1:column
  load(sprintf('chebSigma_%d.mat',i), 'chebSigma')
  cheb_dsdl = diff(chebSigma);
  save(sprintf('cheb_dsdl_%d.mat',i), 'cheb_dsdl')
end

% derivative of sigma with respect to lambda (evaluated)
dsdl = zeros(row,column);
for i=1:column
  load(sprintf('cheb_dsdl_%d.mat',i), 'cheb_dsdl')
  dsdl(:,i) = cheb_dsdl(lambda(:,i));
end
dsdl = reshape(dsdl,[1,M]);

% evaluating sigma at lambda
sl = zeros(row,column);
for i=1:column
  load(sprintf('chebSigma_%d.mat',i), 'chebSigma')
  sl(:,i) = chebSigma(lambda(:,i));
end
sl = reshape(sl,[1,M]);

% repeating minus 1
minusone = repelem(-1,M);

% creating and inverting matrix D
invD = zeros(4,M);
for i=1:M
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

% computing Psi_n

% initializing arrays
F_kT = zeros(row,column);
F_kB = zeros(row,column);
F_kTdl = zeros(row,column);
F_kBdl = zeros(row,column);
F_kTsum = zeros(n,M);
F_kBsum = zeros(n,M);
Psi_nT = zeros(1,M);
Psi_nB = zeros(1,M);

% initializing F_kT and F_kB
for i=1:column
  load(sprintf('chebPsi0BOL_%d.mat',i), 'chebPsi0BOL');
  load(sprintf('chebPsi0TOL_%d.mat',i), 'chebPsi0TOL');
  cheb_F_kTOL = chebPsi0TOL;
  cheb_F_kBOL = chebPsi0BOL;
  save(sprintf('F_kTOL_%d.mat',i), 'cheb_F_kTOL');
  save(sprintf('F_kBOL_%d.mat',i), 'cheb_F_kBOL');
end
term1store = zeros(n,1);
% main loop
for i=1:n
  % computing first term in sumation
  term1 = ((1-sl).^i)/factorial(i);
  % computing F_kT, F_kTdl and F_kBdl
  for j=1:column
    l0 = lambda(1,j);           % begining of domain
    lf = lambda(row,j);         % end of domain

    %loading overlapping chebfuns
    load(sprintf('F_kTOL_%d.mat',j), 'cheb_F_kTOL');
    load(sprintf('F_kBOL_%d.mat',j), 'cheb_F_kBOL');

    % restricting the domains of the overlapping chebfuns
    chebF_kT = restrict(cheb_F_kTOL, [l0,lf]);        % restricted functions
    chebF_kB = restrict(cheb_F_kBOL, [l0,lf]);

    % computing derivative of restricted functions
    chebF_kTdl = diff(chebF_kT);
    chebF_kBdl = diff(chebF_kB);

    % evaluating F_kT at lambda
    F_kT(:,j) = chebF_kT(lambda(:,j));

    F_kTdl(:,j) = chebF_kTdl(lambda(:,j));  % evaluating F_kTdl
    F_kBdl(:,j) = chebF_kBdl(lambda(:,j));  % evaluating F_kBdl
  end

  % reshaping
  F_kT = reshape(F_kT, [1,M]);
  F_kB = reshape(F_kB, [1,M]);
  F_kTdl = reshape(F_kTdl, [1,M]);
  F_kBdl = reshape(F_kBdl, [1,M]);

  % k times top element
  KF_kT = j.*F_kT;

  % computing F_k
  NOinvD_F_k = [F_kTdl(1,:); KF_kT(1,:)+F_kBdl(1,:)];
  F_kT = invD(1,:).*NOinvD_F_k(1,:)+invD(2,:).*NOinvD_F_k(2,:);
  F_kB = invD(3,:).*NOinvD_F_k(1,:)+invD(4,:).*NOinvD_F_k(2,:);

  % storing sumation for each n
  F_kTsum(i,:) = term1.*F_kT;
  F_kBsum(i,:) = term1.*F_kB;

  % updating F_kT and F_kB to F_kTOL and F_kBOL

  % resizing
  F_kT = reshape(F_kT, [row,column]);
  F_kB = reshape(F_kB, [row,column]);
  F_kTdl = reshape(F_kTdl, [row,column]);
  F_kBdl = reshape(F_kBdl, [row,column]);

  % creating new chebfuns with overlaping domains for next iteration
  F_kTOL = [zT;F_kT;zB];
  F_kBOL = [zT;F_kB;zB];

  F_kTOL(repB,column) = F_kTOL(rowOL-domOL,column);
  F_kBOL(repB,column) = F_kBOL(rowOL-domOL,column);
  F_kTOL(repT,1) = F_kTOL(domOL+1,1);
  F_kBOL(repT,1) = F_kBOL(domOL+1,1);

  for j=1:column-1
    F_kTOL(repB,j) = F_kTOL(idB,j+1);   % F_kTOL
    F_kTOL(repT,j+1) = F_kTOL(idT,j);
    F_kBOL(repB,j) = F_kBOL(idB,j+1);   % F_kBOL
    F_kBOL(repT,j+1) = F_kBOL(idT,j);
  end

  % fitting chebfuns to F_kTOL and F_kBOL(as functions of lambda)
  for j=1:column
    l0_OL = lambdaOL(1,j);              % begining of overlaping domain
    lf_OL = lambdaOL(rowOL,j);          % end of overlaping domain

    % new overlaping chebfuns
    cheb_F_kTOL = polyfit(lambdaOL(:,j),F_kTOL(:,j),chebN,domain(l0_OL,lf_OL));
    cheb_F_kBOL = polyfit(lambdaOL(:,j),F_kBOL(:,j),chebN,domain(l0_OL,lf_OL));

    % saving for next iteration
    save(sprintf('F_kTOL_%d.mat',j), 'cheb_F_kTOL');
    save(sprintf('F_kBOL_%d.mat',j), 'cheb_F_kBOL');
  end
end

for i=1:column
  load(sprintf('chebPsi0T_%d.mat',i), 'chebPsi0T');
  load(sprintf('chebPsi0B_%d.mat',i), 'chebPsi0B');
  Psi0T(:,i) = chebPsi0T(lambda(:,i));
  Psi0B(:,i) = chebPsi0B(lambda(:,i));
end

Psi0T = reshape(Psi0T,[1,M]);
Psi0B = reshape(Psi0B,[1,M]);
lambda = reshape(lambda,[1,M]);

Psi_nT = Psi0T+sum(F_kTsum(ChebDiv:n,:));
Psi_nB = Psi0B+sum(F_kBsum(ChebDiv:n,:));

fix = ceil(0.015*M);
Psi_nT(1:fix) = 0;
Psi_nT(M-fix:M) = 0;
Psi_nB(1:fix) = 1;
Psi_nB(M-fix:M) = 1;
Psi_nT_vert = smooth(Psi_nT);
Psi_nB_vert = smooth(Psi_nB);

Psi_nT = Psi_nT_vert.';   % transposing results
Psi_nB = Psi_nB_vert.';


% figure(2)
% plot(lambda, Psi_nT,'LineWidth',2)
% title('Psi_nT')
% axis([0 30 -0.1 0.1])

% figure(3)
% plot(lambda, Psi_nB,'LineWidth',2)
% title('Psi_nB')
% axis([0 30 0.97 1.03])

for i=1:column
  delete(sprintf('cheb_dsdl_%d.mat',i));
  delete(sprintf('chebPsi0B_%d.mat',i));
  delete(sprintf('chebPsi0BOL_%d.mat',i));
  delete(sprintf('chebPsi0T_%d.mat',i));
  delete(sprintf('chebPsi0TOL_%d.mat',i));
  delete(sprintf('chebSigma_%d.mat',i));
  delete(sprintf('chebSigmaOL_%d.mat',i));
  delete(sprintf('F_kBOL_%d.mat',i));
  delete(sprintf('F_kTOL_%d.mat',i));
end


%----------------HODOGRAPH SOLVER------------------%

%-----------INITIALIZING---------%

dLam = lam(2)-lam(1);
dSig = sig(2)-sig(1);

phi = zeros(numLam,numSig);
psi = zeros(numLam,numSig);

courant = dLam/dSig

%-------------BOUNDARY CONDITIONS------------%

phi_bc = Psi_nT; % pulling boundary data
psi_bc = Psi_nB;

phi(:,end) = u;  % implementing boundary conditions at end of domain
psi(:,end) = eta-1;


%psi(:,end) = -0.003*sin(4*lam);        %sin BCs for testing
%phi(:, end) = -0.0000001*sin(40*lam);

%------------INITIAL CONDITIONS---------------%

%phi_ic = 0.00001*sin(sig);

%phi_ic = H1*exp(-c1*(sig - x1).^2) - H2*exp(-c2*(sig - x2).^2);
%psi_ic = H1*exp(-c1*(sig - x1).^2) - H2*exp(-c2*(sig - x2).^2);

%phi(1,:) = phi_ic;
%psi(1,:) = psi_ic;

%--------------MAIN LOOP AND COURANT--------------%

courantList = zeros(numLam, numSig);

for i = 1:numLam-1

   for j = 1:numSig

       %central difference (order O(h)) in space

       if j == numSig

         phi(i+1,j) = phi(i+1, j-1);
         psi(i+1,j) = psi(i+1, j);

       elseif j == 1
         phi(i+1,j) = 0;
         psi(i+1,j) = phi(i,j); %took out sigma
       else
         phi(i+1,j) = ( (-dLam/dSig)*( psi(i,j+1) - psi(i,j-1) ) )+phi(i,j);
         psi(i+1,j) =  ( ( -sig(j)*(dLam/dSig) ) * ( phi(i,j+1) - phi(i,j-1) ) )-(phi(i,j)*dLam)+psi(i,j); %took out sigma
       end

       courantList(i, j) = phi(i+1, j) * courant;

   end
end


%---------------PLOTTING------------------%

figure(1)
mesh(sig, lam, phi)
title('$$\phi$$ in the hodograph plane','interpreter','latex')
xlabel('$$\sigma$$','interpreter','latex')
ylabel('$$\lambda$$','interpreter','latex')
zlabel('$$\varphi$$','interpreter','latex')

figure(2);
mesh(sig, lam, psi);
title('$$\psi$$ in the hodograph plane','interpreter','latex')
xlabel('$$\sigma$$','interpreter','latex')
ylabel('$$\lambda$$','interpreter','latex')
zlabel('$$\psi$$','interpreter','latex')

fprintf('Graphing...\n')
fprintf('Done.\n')

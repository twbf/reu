%----------PDE solver in Hodograph---------%

% This solver uses a finite difference method to solve
% the NSWE in the (\sigma,\lambda) hodograph.

clc; close all; clear all; format Long
warning('off','all')
%load ('strProj.mat','proj')% where boundary data is stored
%load('str.mat','shelf')
addpath('chebfun-master')

%--------------SETUP PARAMTERS------------------%

numSig = 100;
numLam = 4000;

M = 1000; % time steps in runwave.m
N = 1000; % spatial steps in runwave.m


sig = linspace(0,20,numSig); %sigma = 1 is boundary we care about
lam = linspace(0,2,numLam); %lambda, leave at 0-10

Psi_nT = -0.0000001*sin(40*lam);
Psi_nB = -0.005*sin(4*lam);

%-----------INITIALIZING---------%
phi_bc = Psi_nT;
psi_bc = Psi_nB;

dLam = lam(2)-lam(1);
dSig = sig(2)-sig(1);

phi = zeros(numLam,numSig);
psi = zeros(numLam,numSig);

%phi(:, end) = Psi_nT;
psi(:, end) = Psi_nB;


courant = dLam/dSig
Cmin = 0.5;
Cmax = 1;

%-------------BOUNDARY CONDITIONS------------%

%sigVal = 1;
%[~,boundIndex] = min(abs(sig-sigVal)); %locates index for sigma=1

  %phi(:, boundIndex) = phi_bc;
  %psi(:, boundIndex) = psi_bc;

%  phi(:, boundIndex) = zeros(1,numSig);
%  psi(:, boundIndex) = zeros(1,numSig);


%------------INITIAL CONDITIONS---------------%

H1 = 0.006;
H2 = 0;
c1 = 0.44;
c2 = 0;
x1 = 10;
x2 = 0;

phi_ic = 0.00001*sin(sig);
psi_ic = H1*exp(-c1*(sig - x1).^2);

% phi_ic = zeros(1,numSig);
% psi_ic = zeros(1,numSig);

%phi(1,:) = phi_ic;
%phi(2,:) = phi_ic;
psi(1,:) = psi_ic;
psi(2,:) = psi_ic;


fprintf('Please wait...\n');
pause(2)
fprintf('Taking derivatives...\n')

%-----------MAIN LOOP AND COURANT CONDITION------------%

courantList = zeros(numLam, numSig);


for i = 2:numLam-1

   for j = 1:numSig

       %2-point forward difference (order O(h)) both time and space

       if j == numSig

         phi(i+1,j) = phi(i+1, j-1);
         psi(i+1,j) = psi(i+1, j);

       elseif j == 1
         phi(i+1,j) = 0;
         psi(i+1,j) = phi(i+1,j+1); %took out sigma

       else
         %phi(i+1,j) = ( (dLam/dSig) * ( -3*psi(i,j-1) + 4*psi(i,j) - psi(i,j+1) ) ) -3*phi(i-1,j) + 4*phi(i,j);
         %psi(i+1,j) = ( (dLam/dSig*sig(j)) * ( -3*phi(i,j-1) + 4*phi(i,j) - phi(i,j+1)) ) + ( 2*dLam*phi(i,j) ) - 3*psi(i-1,j) + 4*psi(i,j);

         %phi(i+1,j) = ( (dLam/dSig) * ( -3*psi(i,j-1) + 4*psi(i,j) - psi(i,j+1) ) ) + phi(i,j);
         %psi(i+1,j) = ( (dLam/dSig*sig(j)) * ( -3*phi(i,j-1) + 4*phi(i,j) - phi(i,j+1)) ) + ( dLam*phi(i,j) ) + psi(i,j);

         %phi(i+1,j) = (((dLam/dSig)*(3*psi(i,j)-4*psi(i,j+1)+psi(i,j+2)))+phi(i+2,j)+3*phi(i,j))/4;
%psi(i+1,j) = (((dLam/dSig*-sig(i))*(-3*phi(i,j)+4*phi(i,j+1)-phi(i,j+2)))-(2*dLam*phi(i,j))+3*psi(i,j)+psi(i+2,j))/4;
         %phi(i+1,j) = ( (-dLam/dSig)*( psi(i,j+1) - psi(i,j) ) )+phi(i,j);
         %psi(i+1,j) =  ( ( -sig(j)*(dLam/dSig) ) * ( phi(i,j+1) - phi(i,j) ) )-(phi(i,j)*dLam)+psi(i,j); %took out sigma

         phi(i+1,j) = ( (-dLam/dSig)*( psi(i,j+1) - psi(i,j-1) ) )+phi(i,j);
         psi(i+1,j) =  ( ( -sig(j)*(dLam/dSig) ) * ( phi(i,j+1) - phi(i,j-1) ) )-(phi(i,j)*dLam)+psi(i,j); %took out sigma
       end



       courantList(i, j) = phi(i+1, j) * courant;
       %courantList(j) = abs(phi_bc(j))*courant; % courant test


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

figure(3);
mesh(courantList);

fprintf('Graphing...\n')
pause(0.25)
fprintf('Done.\n')

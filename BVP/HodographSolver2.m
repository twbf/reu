%----------PDE solver in Hodograph---------%

% This solver uses a finite difference method to solve 
% the NSWE in the (\sigma,\lambda) hodograph.

clc; close all; clear all; format Long
warning('off','all')
load ('strProj.mat','proj')% where boundary data is stored
load('str.mat','shelf')
addpath('chebfun-master')

%--------------SETUP PARAMTERS------------------%

M = shelf.M; % time steps in runwave.m
N = shelf.N; % spatial steps in runwave.m

Psi_nT = proj.phi;
Psi_nB = proj.psi;

numSig = length(Psi_nB);
numLam = length(Psi_nT);

sig = linspace(0,8,numSig); %sigma = 1 is boundary we care about
lam = linspace(0,0.006,numLam); %lambda, leave at 0-10

%-----------INITIALIZING---------%
phi_bc = Psi_nT;
psi_bc = Psi_nB;

dLam = lam(2)-lam(1);
dSig = sig(2)-sig(1);

phi = zeros(numLam,numSig);
psi = zeros(numLam,numSig);

courant = dLam/dSig
Cmin = 0.5; 
Cmax = 1;

%-------------BOUNDARY CONDITIONS------------%

sigVal = 1;
[~,boundIndex] = min(abs(sig-sigVal)); %locates index for sigma=1

  phi(:, boundIndex) = phi_bc;
  psi(:, boundIndex) = psi_bc;

%  phi(:, boundIndex) = zeros(1,numSig);
%  psi(:, boundIndex) = zeros(1,numSig);


%------------INITIAL CONDITIONS---------------%

H1 = 0.6;
H2 = 0;
c1 = 0.4444;
c2 = 0;
x1 = 4.1209;
x2 = 0;

phi_ic = H1*exp(-c1*(sig - x1).^2);
psi_ic = H1*exp(-c1*(sig - x1).^2);

% phi_ic = zeros(1,numSig);
% psi_ic = zeros(1,numSig);

phi(1,:) = phi_ic;
psi(1,:) = psi_ic;
 

fprintf('Please wait...\n');
pause(2)
fprintf('Taking derivatives...\n')

%-----------MAIN LOOP AND COURANT CONDITION------------%

courantList(M) = abs(phi_bc(end))*courant;


for i = 1:numLam-2

   for j = 1:numLam-2

       %2-point forward difference (order O(h)) both time and space

       phi(i+1,j) = ((-dLam/dSig)*(psi(i,j+1)-psi(i,j)))+phi(i,j);
       psi(i+1,j) =  ((-sig(i)*(dLam/dSig))*(phi(i,j+1)-phi(i,j)))-(phi(i,j)*dLam)+psi(i,j);
    
    courantList(j) = abs(phi_bc(j))*courant; % courant test 
    
    
   end
end


%---------------PLOTTING------------------%

figure(1)
mesh(sig, lam, phi)
title('$$\phi$$ in the hodograph plane','interpreter','latex')
xlabel('$$\sigma$$','interpreter','latex')
ylabel('$$\lambda$$','interpreter','latex')
zlabel('$$\varphi$$','interpreter','latex')


figure(2)
mesh(sig, lam, psi)
title('$$\psi$$ in the hodograph plane','interpreter','latex')
xlabel('$$\sigma$$','interpreter','latex')
ylabel('$$\lambda$$','interpreter','latex')
zlabel('$$\psi$$','interpreter','latex')

fprintf('Graphing...\n')
pause(0.25)
fprintf('Done.\n')





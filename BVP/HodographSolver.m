%----------PDE solver in Hodograph---------%

% This solver uses ??? Method to solve the NSWE in the (\sigma,\lambda)
% hodograph.

clc; close all; clear all; format Long
warning('off','all')
%load ('Psi_n.mat','Psi_nB','Psi_nT')
addpath('chebfun-master')

numSig = 100;
numLam = 10000;

sig = linspace(0,1,numSig); %sigma
lam = linspace(0,3,numLam); %lambda

Psi_nB = 0.001*cos(lam);
Psi_nT = 0.0000001*sin(lam);

%--------Defining the length of our vectors------%


%----Loading Boundary Data and creating vectors---%
phi_bc = Psi_nT;
psi_bc = Psi_nB;

dSig = sig(2)-sig(1);
dLam = lam(2)-lam(1);

phi = zeros(numLam,numSig);
psi = zeros(numLam,numSig);


%---------Courant test for stability----------%
courant = zeros(1,numLam);
deltaX = dLam;  % what's on the x-axis
deltaT = dSig;   % what's on the y-axis

for i=1:numLam
  courant(i) = (abs(phi_bc(i))*deltaT)/deltaX;
end

courant = sort(courant);
courant(end)
if courant >=0.5 | courant <=0.1
  disp('This is potentially unstable.')
end

%------------Initializing vectors----------%

sigVal = 1;
[~,boundIndex] = min(abs(sig-sigVal)); %locates index for sigma=1

phi_ic = zeros(1,numSig);
psi_ic = zeros(1,numSig);

phi(1,:) = phi_ic;
psi(1,:) = psi_ic;

phi(:, boundIndex) = phi_bc;
psi(:, boundIndex) = psi_bc;


for i = 2:numLam-1

   for j = 2:numSig-1

       % Two-point central difference in space
       % Two-point central difference in time 0(h^2)

    phi(i+1,j) = ((-dLam/dSig)*(psi(i,j+1)-psi(i,j-1)))+phi(i-1,j);
    psi(i+1,j) = -((dLam/dSig)*(phi(i,j+1)-phi(i,j-1)))-(phi(i,j)*2*dLam)+psi(i-1,j);

    % note: should be a c-squared sigma term in psi

   end
end

%mesh(x,y,z) = (lambda,sigma,function)

figure(1)
mesh(sig, lam, phi)
title('$$\phi$$ in the hodograph plane','interpreter','latex')
ylabel('$$\lambda$$','interpreter','latex')
xlabel('$$\sigma$$','interpreter','latex')
zlabel('$$\varphi$$','interpreter','latex')

figure(2)
mesh(sig, lam, psi)
title('$$\psi$$ in the hodograph plane','interpreter','latex')
ylabel('$$\lambda$$','interpreter','latex')
xlabel('$$\sigma$$','interpreter','latex')
zlabel('$$\psi$$','interpreter','latex')

%----------PDE solver in Hodograph---------%

% This solver uses a finite difference method to solve
% the NSWE in the (\sigma,\lambda) hodograph.

clc; close all; clear all; format Long
warning('off','all')
load('strProj.mat','proj')
addpath('chebfun-master')



%--------------SETUP PARAMTERS------------------%

numSig = 5;
numLam = 10000;

sig = linspace(0,1,numSig); %sigma = 1 is boundary we care about
lam = linspace(0,10,numLam); %lambda, leave at 0-10

%-----------INITIALIZING---------%

dLam = lam(2)-lam(1);
dSig = sig(2)-sig(1);

phi = zeros(numLam,numSig);
psi = zeros(numLam,numSig);

courant = dLam/dSig
Cmin = 0.5;
Cmax = 1;

%-------------BOUNDARY CONDITIONS------------%

%  phi_bc = proj.phi;
%  psi_bc = proj.psi;
 
%psi(:,end) = -0.003*sin(4*lam);

%phi(:, end) = -0.0000001*sin(40*lam);

psi_bc = -0.003*sin(4*lam);


%------------INITIAL CONDITIONS---------------%

H1 = 0.006;
c1 = 0.44;
x1 = 5;

phi_ic = 0.00001*sin(sig);
psi_ic = H1*exp(-c1*(sig - x1).^2);

%phi(1,:) = phi_ic;
psi(1,:) = psi_ic;

fprintf('Please wait...\n');
fprintf('Taking derivatives...\n')

%-----------MAIN LOOP AND COURANT CONDITION------------%

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

% figure(3);
% mesh(courantList);

fprintf('Graphing...\n')
fprintf('Done.\n')




%----------PDE solver in Hodograph---------%

% This solver uses a finite difference method to solve
% the NSWE in the (\sigma,\lambda) hodograph.


function [phi,psi] = HodoSolve(Psi)

  global g td
  global t0 Tf t_res x_res

  %--------------SETUP PARAMTERS------------------%

  numSig = 15;

  sig = linspace(0,1,numSig); %sigma = 1 is boundary we care about
  lam = linspace(t0,Tf,t_res); %lambda, leave at 0-10

  %-----------INITIALIZING---------%

  dLam = lam(2)-lam(1);
  dSig = sig(2)-sig(1);
  courant = dLam/dSig;

  phi = zeros(t_res,numSig);
  psi = zeros(t_res,numSig);

  %-------------BOUNDARY CONDITIONS------------%

  phi(:, end) = Psi(1,:);
  psi(:, end) = Psi(2,:);
  %psi_bc = -0.003*sin(4*lam);

  %------------INITIAL CONDITIONS---------------%

  %phi(1,:) = 0.00001*sin(sig);
  %psi(1,:) = H1*exp(-c1*(sig - x1).^2);



  fprintf('Please wait...\n');
  fprintf('Taking derivatives...\n')

%-----------MAIN LOOP AND COURANT CONDITION------------%

courantList = zeros(t_res, numSig);

for i = 1:t_res-1

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

figure(5)
mesh(sig, lam, phi)
title('Solution, $$\phi$$ in the hodograph plane','interpreter','latex')
xlabel('$$\sigma$$','interpreter','latex')
ylabel('$$\lambda$$','interpreter','latex')
zlabel('$$\varphi$$','interpreter','latex')


figure(6);
mesh(sig, lam, psi);
title('Solution, $$\psi$$ in the hodograph plane','interpreter','latex')
xlabel('$$\sigma$$','interpreter','latex')
ylabel('$$\lambda$$','interpreter','latex')
zlabel('$$\psi$$','interpreter','latex')

% figure(3);
% mesh(courantList);

fprintf('Graphing...\n')
fprintf('Done.\n')


end

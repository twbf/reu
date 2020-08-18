%----------PDE solver in Hodograph---------%

% This solver uses a finite difference method to solve
% the NSWE in the (\sigma,\lambda) hodograph.


function [eta,u] = HodoSolve(Psi)

  global g td
  global t0 Tf t_res x_res numSig

  %--------------SETUP PARAMTERS------------------%

  sig = linspace(0,1,numSig); %sigma = 1 is boundary we care about
  lam = linspace(t0,Tf*sqrt(g),t_res); %lambda, leave at 0-10

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

for i = 2:t_res-1

   for j = 1:numSig

       %central difference (order O(h)) in space

       if j == numSig

         %phi(i+1,j) = phi(i+1, j);
         %psi(i+1,j) = psi(i+1, j);

       elseif j ==1

        %psi(i+1,j) = -(phi(i,j+1)-phi(i,j))*( (dLam)^2/(dSig*2) ) - psi(i-1,j) +2*psi(i,j);
        %phi(i+1,j) = -(psi(i,j+1)-psi(i,j))*( (dLam)^2/(dSig*2) ) - phi(i-1,j) +2*phi(i,j);
         %phi(i+1,j) = 0;
         %psi(i+1,j) = psi(i,j-1); %took out sigma
         %psi(i+1,j) =  ( ( -3.5*sig(j+1)*(dLam/dSig) ) * ( phi(i,j+1) - phi(i,j) ) )-(phi(i,j)*2*dLam)+psi(i-1,j);
       else

         phi(i+1,j) = ( (-dLam/dSig)*( psi(i,j+1) - psi(i,j-1) ) )+phi(i-1,j);

         psi(i+1,j) =  ( ( -sig(j)*(dLam/dSig) ) * ( phi(i,j+1) - phi(i,j-1) ) )-(phi(i,j)*2*dLam)+psi(i-1,j);
         
       end

       courantList(i, j) = phi(i+1, j) * courant;

   end
   %solving with simple 2 point FD
   %phi(i+1,1) = ( -(dLam/dSig) * ( psi(i+1,1+1) - psi(i,1) ) + phi(i,1) )/( 1 + (dLam)^2 / dSig );
   %psi(i+1,1) = -dLam * phi(i+1,1) + psi(i,1);

   %2 point central differnce
   %phi(i+1,1) = ( -(2*dLam/dSig) * ( psi(i+1,1+1) - psi(i,1) ) + phi(i-1,1) )/( 1 + (2*dLam)^2 / dSig );
   %psi(i+1,1) = -dLam *2* phi(i+1,1) + psi(i-1,1);

   %2 point without solving
   psi(i+1,1) = -dLam*phi(i,1) + psi(i,1);
   phi(i+1,1) = ( (-dLam/dSig)*( psi(i,1+1) - psi(i,1) ) )+phi(i,1);

   %3 point backwards difference without solving
   %phi(i+1,1) = ( ( (-2*dLam/dSig)*( psi(i,1+1) - psi(i,1) ) )-phi(i-1,1) + 4*phi(i,1) )/3;
   %psi(i+1,1) = (-2*dLam*phi(i+1,1) - psi(i-1,1) + 4*psi(i,1))/3;

   %3 point backwards differnce
   %phi(i+1,1) = ( -(2*dLam/dSig) * ( psi(i+1,1+1) - psi(i,1) ) + phi(i-1,1) )/( 1 + (2*dLam)^2 / dSig );
   %psi(i+1,1) = -dLam *2* phi(i+1,1) + psi(i-1,1);

   %phi(i+1,1) = phi(i,1)/(1+(dLam)^2/dSig);

   %secound derivative heat
   %psi(i+1,1) = ( psi(i+1,1+1) - dSig/(dLam)^2*( psi(i-1,1) - 2*psi(i,1) ) )/( 1 + dSig/(dLam)^2 );
   %phi(i+1,1) = ( phi(i+1,1+1) - dSig/(dLam)^2*( phi(i-1,1) - 2*phi(i,1) ) )/( 1 + dSig/(dLam)^2 );


   %guess
   %psi(i+1,1) = -(phi(i+1,1+1)-phi(i,1))*( (dLam)^2/(dSig*2) ) - psi(i-1,1) +2*psi(i,1);

   %phi(i+1,1) = (psi(i+1,1+1)-psi(i+1,1))*( (dLam)^2/(dSig*2) ) - phi(i-1,1) +2*phi(i,1);
   %psi(i+1,1) = -(phi(i+1,1+1)-phi(i+1,1))*( (dLam)^2/(dSig*2) ) - psi(i-1,1) +2*psi(i,1);
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

  %CG transform

  disp('    backwards CG transform and demensionalization... ');

  eta = zeros(t_res, numSig);
  u = zeros(t_res, numSig);
  xx = zeros(t_res, numSig);
  tt = zeros(t_res, numSig);

  for i=1:t_res
    %CG transform
    u(i,:) = phi(i,:);
    eta(i,:) = psi(i,:) - u(i,:).^2/2;
    xx(i,:) = sig - eta(i,:);
    tt(i,:) = u(i,:) + lam(i);

    %deminsionalizing
    u(i,:) = u(i,:)*sqrt(g*td);
    eta(i,:) = eta(i,:);
    tt(i,:) = tt(i,:)/sqrt(td*g);
  end

  %to display eta and u
    %figure(5);
    %mesh(tt,xx,eta);

    %figure(5);
    %mesh(tt,xx,u);

  disp('    scattered interpolation of eta and u... ');

  s_tt = reshape(tt, [t_res*numSig, 1]);
  s_xx = reshape(xx, [t_res*numSig, 1]);
  s_eta = reshape(eta, [t_res*numSig, 1]);
  s_u = reshape(u, [t_res*numSig, 1]);

  eta = scatteredInterpolant(s_xx, s_tt, s_eta);
  u = scatteredInterpolant(s_xx, s_tt, s_u);

end

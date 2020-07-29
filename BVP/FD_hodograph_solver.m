%number of s and lambda points
num_s = 1000;
num_la = 100;

%range over which we solve
s = linspace(0, 10, num_s);
la = linspace(0, 1, num_la);

ds = s(2) - s(1);
dla = la(2) - la(1);

%initial conditions are zero right now
%initial condition
phi_0 = s*0;
psi_0 = s*0;

%boundary
phi_b = la*0;
psi_b = la*0 + -0.000000001;

phi = zeros(num_la, num_s);
psi = zeros(num_la, num_s);

%setting IC
phi(1,:) = phi_0;
psi(1,:) = psi_0;

%setting BC
phi(:,end) = phi_b';
psi(:,end) = psi_b';

for i = 1:(num_la-1)

  for j = 1:(num_s-1)

    phi(i+1, j) = -dla*(psi(i, j+1) - psi(i, j))/ds + phi(i, j);
    psi(i+1, j) = -dla*(phi(i,j) + (phi(i,j+1) - phi(i,j))/ds)  + psi(i,j);

  end

end

figure(1);
mesh(phi);

figure(2);
mesh(psi);

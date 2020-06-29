% Goal: For a given phi(s, lambda) and psi(s, lambda) compute eta, u, x, t


load('psi_phi_old.mat')

%Parameters -----
size = 100;  %resolution
lambda_list = linspace(0, 9, size); %lambda space
s = linspace(0.1, 9, size); %s space   note: at s = 0 regularization is needed - more on this later

g = sqrt(9.81);

x = zeros(size, size);
t = zeros(size, size);
h = zeros(size, size);
u = zeros(size, size);

for i=1:size
    lambda = lambda_list(i); %for a single lambda
    
    %CG transform from (s, lambda) to (x,t)
    u(i,:) = psi(s,lambda);
    h(i,:) = phi(s,lambda)-u(i,:).^2/2;
    x(i,:) = s - h(i,:);
    t(i,:) = (u(i,:) + lambda)/g;  % deminsionalizing
end

% scatter() and scatteredInterpolant() takes lists not matricies so x and t have to be reshaped
xx = reshape(x, [size*size,1]);
tt = reshape(t, [size*size,1]);
hh = reshape(h, [size*size,1]);
uu = reshape(u, [size*size,1]);

figure(1);
scatter(xx,tt);
xlabel('$x$', 'interpreter', 'LaTeX', 'fontsize', 15);
ylabel('$t$', 'interpreter', 'LaTeX', 'fontsize', 15);
title('Equally Spaced Grid in $(s, \lambda)$ Transformed to $(x, t)$ ', 'interpreter', 'LaTeX', 'fontsize', 20);
export_fig('CG_coordinate_transform', '-m2', '-a4', '-painters');

figure(2);
mesh(x,t,h);
xlabel('$x$', 'interpreter', 'LaTeX', 'fontsize', 15);
ylabel('$t$', 'interpreter', 'LaTeX', 'fontsize', 15);
zlabel('$\eta$', 'interpreter', 'LaTeX', 'fontsize', 15);
title(' CG transformed $\eta$', 'interpreter', 'LaTeX', 'fontsize', 20);
export_fig('CG_eta', '-m2', '-a4', '-painters');

figure(3);
mesh(x,t,u);
xlabel('$x$', 'interpreter', 'LaTeX', 'fontsize', 15);
ylabel('$t$', 'interpreter', 'LaTeX', 'fontsize', 15);
zlabel('$u$', 'interpreter', 'LaTeX', 'fontsize', 15);
title('CG transformed $u$', 'interpreter', 'LaTeX', 'fontsize', 20);
export_fig('CG_u', '-m2', '-a4', '-painters');


%interpolating mesh in order to get a smooth polynomial curve
analytic_eta = scatteredInterpolant(xx,tt,hh);
save('analytic_eta', 'analytic_eta');


analytic_u = scatteredInterpolant(xx,tt,u);
save('analytic_u', 'analytic_u');


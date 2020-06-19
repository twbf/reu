clc
clear
close all
format longE

addpath('export_fig/');

%%% Physical parameter:
x0 = 5;
K  = 4.0;		% upper bound in k domain [0, K]
L  = 20.0;		% upper bound in x domain [0, L]
Ls = 10.0;		% upper bound for s parameter
La = 10.0;		% upper bound for lambda parameter

%%% SWE Parameters

m = Inf
beta = 1

disp('eta....')

x   = chebfun('x', [0 L]);
eta = chebfun(@(x) exp(-(x - x0).^2), [0 L]);
eta_prime = chebfun(@(x) -2*(x-x0)*exp(-(x - x0).^2), [0 L]);

disp('done')

 figure(1);
 plot(eta, '-', 'LineWidth', 2.0), grid off
 xlabel('$x$', 'interpreter', 'LaTeX', 'fontsize', 12);
 ylabel('$\eta(x)$', 'interpreter', 'LaTeX', 'fontsize', 12);
 title('Initial free surface elevation');
 set(gcf, 'color', 'w');
 export_fig('Eta.png', '-m2', '-a4', '-painters');


disp('j0....')

j0 = chebfun(@(kx) besselj(0.0, kx), [0.0 max([2.0*K*max(sqrt(x + eta)) 2.0*K*sqrt(Ls)])]);

j1 = chebfun(@(kx) besselj(1.0, kx), [0.0 max([2.0*K*max(sqrt(x + eta)) 2.0*K*sqrt(Ls)])]);

disp('a....')
a  = chebfun(@(k) 2*k*sum(eta(x).*j0(2.0*k*sqrt(x + eta(x))).*(1 + eta_prime(x))), [0 K]);

disp('done')

 figure(2);
 plot(a, '-', 'LineWidth', 2.0), grid off
 xlabel('$k$', 'interpreter', 'LaTeX', 'fontsize', 12);
 ylabel('$a(k)$', 'interpreter', 'LaTeX', 'fontsize', 12);
 title('Coefficient $a(k)$', 'interpreter', 'LaTeX', 'fontsize', 13);
 set(gcf, 'color', 'w');
 export_fig('Eta1.png', '-m2', '-a4', '-painters');

disp('phi....')

Cos = chebfun(@(lk) cos(lk), [0 La*K], 'vectorize');
k   = chebfun('x', [0 K]);
phi   = chebfun2(@(s,la) sum(a(k)*Cos(la*k)*j0(2.0*k*sqrt(s))), [0 Ls 0 La], 'vectorize');


 figure(3);
 plot(phi);
 xlabel('$s$', 'interpreter', 'LaTeX', 'fontsize', 12);
 ylabel('$\lambda$', 'interpreter', 'LaTeX', 'fontsize', 12);
 view([0 90]); colorbar;
 title('Phi Two-parameters integral $f(s,\lambda)$', 'interpreter', 'LaTeX', 'fontsize', 12);
 set(gcf, 'color', 'w');
 export_fig('phi.png', '-m2', '-a4', '-painters');

 disp('done')
disp(phi(1,2))

 disp('psi...')

Sin = chebfun(@(lk) sin(lk), [0 La*K], 'vectorize');
psi = chebfun2(@(s,la) s^(-1/2)*sum(a(k)*Sin(la*k)*j1(2.0*k*sqrt(s))), [0.01 Ls 0.01 La], 'vectorize');

disp('done')

 figure(4);
 plot(psi);
 xlabel('$s$', 'interpreter', 'LaTeX', 'fontsize', 12);
 ylabel('$\lambda$', 'interpreter', 'LaTeX', 'fontsize', 12);
 view([0 90]); colorbar;
 title('Psi Two-parameters integral $f(s,\lambda)$', 'interpreter', 'LaTeX', 'fontsize', 12);
 set(gcf, 'color', 'w');
 export_fig('psi.png', '-m2', '-a4', '-painters');

save('psi_phi')



%u = psi(s,lamda)

%eta = phi(s,lambda)-u^2/2

%x = s-eta

%t = u + lamda

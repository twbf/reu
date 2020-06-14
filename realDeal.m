clc
clear
close all
format longE

addpath('export_fig/');

%%% Physical parameter:
x0 = 3.5;
K  = 20.0;		% upper bound in k domain [0, K]
L  = 20.0;		% upper bound in x domain [0, L]
Ls = 10.0;		% upper bound for s parameter
La = 10.0;		% upper bound for lambda parameter

x   = chebfun('x', [0 L]);
eta = chebfun(@(x) exp(-(x - x0).^2), [0 L]);

% figure(1);
% plot(eta, '-', 'LineWidth', 2.0), grid off
% xlabel('$x$', 'interpreter', 'LaTeX', 'fontsize', 12);
% ylabel('$\eta(x)$', 'interpreter', 'LaTeX', 'fontsize', 12);
% title('Initial free surface elevation');
% set(gcf, 'color', 'w');
% export_fig('shots/Eta.png', '-m2', '-a4', '-painters');

j0 = chebfun(@(kx) besselj(0.0, kx), [0.0 max([2.0*K*max(sqrt(x + eta).*(1 + eta)) 2.0*K*sqrt(Ls)])]);
a  = chebfun(@(k) sum(eta(x).*j0(2.0*k*sqrt(x + eta(x)).*(1 + eta(x)))), [0 K]);

% figure(2);
% plot(a, '-', 'LineWidth', 2.0), grid off
% xlabel('$k$', 'interpreter', 'LaTeX', 'fontsize', 12);
% ylabel('$a(k)$', 'interpreter', 'LaTeX', 'fontsize', 12);
% title('Coefficient $a(k)$', 'interpreter', 'LaTeX', 'fontsize', 13);
% set(gcf, 'color', 'w');
% export_fig('shots/Eta.png', '-m2', '-a4', '-painters');

Cos = chebfun(@(lk) cos(lk), [0 La*K], 'vectorize');
k   = chebfun('x', [0 K]);
F   = chebfun2(@(s,la) sum(a*Cos(la*k)*j0(2.0*k*sqrt(s))), [0.05 Ls 0.05 La], 'vectorize');

% figure(3);
% plot(F);
% xlabel('$s$', 'interpreter', 'LaTeX', 'fontsize', 12);
% ylabel('$\lambda$', 'interpreter', 'LaTeX', 'fontsize', 12);
% view([0 90]); colorbar;
% title('Two-parameters integral $f(s,\lambda)$', 'interpreter', 'LaTeX', 'fontsize', 12);
% set(gcf, 'color', 'w');
% export_fig('shots/MainIntegral.png', '-m2', '-a4', '-painters');

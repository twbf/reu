clear
close all
format longE

addpath('export_fig/');

%%% Example 1:

Kmax = 10.0;	% the range of k variable is [0, Kmax]
L = 10.0; 		% upper integration limit in Hankel transform

a  = 1.0;		% parameter
f  = chebfun(@(x) exp(-0.5*a^2*x^2), [0 L]);
j0 = chebfun(@(kx) besselj(0,kx), [0 L*Kmax]);

% Forward transform:
x = chebfun('x', [0 L]);
F = chebfun(@(k) sum(f.*j0(k.*x).*x), [0 Kmax]);

% Backward transform:
k = chebfun('x', [0 Kmax]);
Finv = chebfun(@(x) sum(F(k).*j0(k.*x).*k), [0 L]);

% Exact transform:
Fexa = @(k) 1.0/(a*a)*exp(-0.5*k.^2/(a*a));
% source: https://en.wikipedia.org/wiki/Hankel_transform

% Forward:

Ks = linspace(0.0, Kmax, 100);
Fe = Fexa(Ks);
Fn = F(Ks);

figure(1);
subplot(2,1,1);
plot(Ks, Fe, 'r-', 'LineWidth', 2.5), grid off, hold on
plot(Ks, Fn, 'b--', 'LineWidth', 2.0), hold off
ll = legend('\ Exact transform', '\ Numerical');
set(ll, 'box', 'off');
set(ll, 'Interpreter', 'LaTeX');
set(ll, 'FontSize', 12);
xlabel('$k$', 'Interpreter', 'LaTeX', 'FontSize', 12);
ylabel('$F_0(k)$', 'Interpreter', 'LaTeX', 'FontSize', 12);
title('Forward Hankel transform: $\exp(-\frac{1}{2}a^2 x^2)$', 'Interpreter', 'LaTeX');

subplot(2,1,2);
plot(Ks, Fe-Fn, 'r-', 'LineWidth', 1.5), grid off, hold off
xlabel('$k$', 'Interpreter', 'LaTeX', 'FontSize', 12);
ylabel('$F_{\rm exact} - F_{\rm num}$', 'Interpreter', 'LaTeX', 'FontSize', 12);
set(gcf, 'color', 'w');

export_fig('shots/HankelTest1.png', '-m2', '-a4', '-painters');

% Backward:

X  = linspace(0.0, L, 100);
fX = f(X);
Fi = Finv(X);

figure(2);
subplot(2,1,1);
plot(X, fX, 'r-', 'LineWidth', 2.5), grid off, hold on
plot(X, Fi, 'b--', 'LineWidth', 2.0), hold off
ll = legend('\ Exact transform', '\ Numerical');
set(ll, 'box', 'off');
set(ll, 'Interpreter', 'LaTeX');
set(ll, 'FontSize', 12);
xlabel('$r$', 'Interpreter', 'LaTeX', 'FontSize', 12);
ylabel('$f(r)$', 'Interpreter', 'LaTeX', 'FontSize', 12);
title('Backward Hankel transform: $\exp(-\frac{1}{2}a^2 x^2)$', 'Interpreter', 'LaTeX');

subplot(2,1,2);
plot(X, fX-Fi, 'r-', 'LineWidth', 1.5), grid off, hold off
xlabel('$r$', 'Interpreter', 'LaTeX', 'FontSize', 12);
ylabel('$F_{\rm exact} - F_{\rm num}$', 'Interpreter', 'LaTeX', 'FontSize', 12);
set(gcf, 'color', 'w');

export_fig('shots/HankelBackTest1.png', '-m2', '-a4', '-painters');

clear

%%% Example 2:

Kmax = 100.0;	% the range of k variable is [0, Kmax]
L = 20.0; 		% upper integration limit in Hankel transform

a  = 2.0;		% parameter
f  = chebfun(@(x) exp(-a*x), [0 L]);
j0 = chebfun(@(kx) besselj(0,kx), [0 L*Kmax]);

% Forward transform:
x = chebfun('x', [0 L]);
F = chebfun(@(k) sum(f.*j0(k.*x).*x), [0 Kmax]);

% Backward transform:
k = chebfun('x', [0 Kmax]);
Finv = chebfun(@(x) sum(F(k).*j0(k.*x).*k), [0 L]);

% Exact transform:
Fexa = @(k) a./((a^2 + k.^2)).^(1.5);

Ks = linspace(0.0, Kmax, 1000);
Fe = Fexa(Ks);
Fn = F(Ks);

% Forward:

figure(3);
subplot(2,1,1);
plot(Ks, Fe, 'r-', 'LineWidth', 2.5), grid off, hold on
plot(Ks, Fn, 'b--', 'LineWidth', 2.0), hold off
xlim([0 Kmax/2]);
ll = legend('\ Exact transform', '\ Numerical');
set(ll, 'box', 'off');
set(ll, 'Interpreter', 'LaTeX');
set(ll, 'FontSize', 12);
xlabel('$k$', 'Interpreter', 'LaTeX', 'FontSize', 12);
ylabel('$F_0(k)$', 'Interpreter', 'LaTeX', 'FontSize', 12);
title('Forward Hankel transform: $\exp(-a\cdot x)$', 'Interpreter', 'LaTeX');

subplot(2,1,2);
plot(Ks, Fe-Fn, 'r-', 'LineWidth', 1.5), grid off, hold off
xlim([0 Kmax/2]);
xlabel('$k$', 'Interpreter', 'LaTeX', 'FontSize', 12);
ylabel('$F_{\rm exact} - F_{\rm num}$', 'Interpreter', 'LaTeX', 'FontSize', 12);
set(gcf, 'color', 'w');

export_fig('shots/HankelTest2.png', '-m2', '-a4', '-painters');

% Backward:
X  = linspace(0.0, L, 1000);
fX = f(X);
Fi = Finv(X);

figure(4);
subplot(2,1,1);
plot(X, fX, 'r-', 'LineWidth', 2.5), grid off, hold on
plot(X, Fi, 'b--', 'LineWidth', 2.0), hold off
ll = legend('\ Exact transform', '\ Numerical');
set(ll, 'box', 'off');
set(ll, 'Interpreter', 'LaTeX');
set(ll, 'FontSize', 12);
xlim([0 L/2]);
xlabel('$r$', 'Interpreter', 'LaTeX', 'FontSize', 12);
ylabel('$f(r)$', 'Interpreter', 'LaTeX', 'FontSize', 12);
title('Backward Hankel transform: $\exp(-a\cdot x)$', 'Interpreter', 'LaTeX');

subplot(2,1,2);
plot(X, fX-Fi, 'r-', 'LineWidth', 1.5), grid off, hold off
xlim([0 L/2]);
xlabel('$r$', 'Interpreter', 'LaTeX', 'FontSize', 12);
ylabel('$F_{\rm exact} - F_{\rm num}$', 'Interpreter', 'LaTeX', 'FontSize', 12);
set(gcf, 'color', 'w');

export_fig('shots/HankelBackTest2.png', '-m2', '-a4', '-painters');
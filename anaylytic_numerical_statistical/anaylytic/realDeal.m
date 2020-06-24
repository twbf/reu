clc
clear
close all
format longE

addpath('export_fig/');

global eta eta_prime u u_prime beta

%%% Physical parameter:
x0 = 5;
K  = 4.0;		% upper bound in k domain [0, K]
L  = 20.0;		% upper bound in x domain [0, L]
Ls = 10.0;		% upper bound for s parameter
La = 10.0;		% upper bound for lambda parameter

%%% SWE Parameters

m = Inf;
beta = 1;

disp('eta, eta_prime, u....')

H1 = 1.006;
H2 = 0.018;
c1 = 0.4444;
c2 = 4.0;
x1 = 4.1209;
x2 = 1.6384;

x   = chebfun('x', [0 L]);
eta = chebfun(@(x) H1*exp(-c1*(x - x1).^2) - H2*exp(-c2*(x - x2).^2), [0 L]);
eta_prime = chebfun(@(x) -2*H1*c1*(x-x1)*exp(-c1*(x - x1).^2) + -2*H2*c2*(x-x2)*exp(-c2*(x - x2).^2), [0 L]);

u   = chebfun(@(x) exp(-(x-x0)^.2), [0 L]);
u_prime = chebfun(@(x) -2*(x-x0)*exp(-(x-x0)^.2), [0 L]);

%Plotting gamma
%     t = linspace(0,10,100);
%     fir = t+eta(t);
%     sec = -u(t);
%     scatter(fir,sec)

%data projection to lambda = 0
%I can make this faster by setting variables


 %figure(1);
 %plot(eta)
 %xlabel('$x$', 'interpreter', 'LaTeX', 'fontsize', 12);
 %ylabel('$\eta(x)$', 'interpreter', 'LaTeX', 'fontsize', 12);
 %title('Initial free surface elevation');
 %set(gcf, 'color', 'w');
 %export_fig('Eta.png', '-m2', '-a4', '-painters');


disp('j0...')
j0 = chebfun(@(kx) besselj(0.0, kx), [0.0 max([2.0*K*max(sqrt(x + eta)) 2.0*K*sqrt(Ls)])]);

disp('j1...')
j1 = chebfun(@(kx) besselj(1.0, kx), [0.0 max([2.0*K*max(sqrt(x + eta)) 2.0*K*sqrt(Ls)])]);

disp('cos...')
Cos = chebfun(@(lk) cos(lk), [0 La*K], 'vectorize');

disp('sin...')
Sin = chebfun(@(lk) sin(lk), [0 La*K], 'vectorize');

disp('p...')
p = chebfun('x', [0 K]);


s = @(x) x + eta(x);
A = @(x) [0 1; beta^2*s(x) 0]; %beta^2*s
B = [0 0; 1 0];
    
D = @(x) eta_prime(x)*eye(2) + u_prime(x)*A(x);
    
phi0 = @(x) [u(x) ; eta(x)+(u(x).^2)/2];
phi0_prime = @(x) [u_prime(x); eta_prime(x)+2*u(x)*u_prime(x)];
    
proj = @(x) phi0(x) + u(x)*(u_prime(x)*inv(D(x))*B*phi0(x) - B*phi0(x) -A(x)*inv(D(x))*phi0_prime(x));

disp("Proj1")
test = linspace(0,1,20);

disp(size(test))

for i=1:size(test,2)
    out(:,:,i) = proj(i);
end

disp(out)


q = integral(proj,0,1);
disp(p)


disp('a...')
a  = chebfun(@(k) 2*k*sum( [0 1]*proj(p)*j0(2*k*sqrt(p)) ), [0 K]);

disp('b...')
b  = chebfun(@(k) -2*beta*k*sum( [1 0]*proj1(p)*p^(1/2)*j1(2*k*sqrt(p)) ), [0 K]);



%if initial velocity = 0

% disp('a....')
% a  = chebfun(@(k) 2*k*sum(eta(x).*j0(2.0*k*sqrt(x + eta(x))).*(1 + eta_prime(x))), [0 K]);
% 
% disp('done')
% 
%  figure(2);
%  plot(a, '-', 'LineWidth', 2.0), grid off
%  xlabel('$k$', 'interpreter', 'LaTeX', 'fontsize', 12);
%  ylabel('$a(k)$', 'interpreter', 'LaTeX', 'fontsize', 12);
%  title('Coefficient $a(k)$', 'interpreter', 'LaTeX', 'fontsize', 13);
%  set(gcf, 'color', 'w');
%  export_fig('Eta1.png', '-m2', '-a4', '-painters');
% 
% disp('phi....')
% 
% Cos = chebfun(@(lk) cos(lk), [0 La*K], 'vectorize');
% k   = chebfun('x', [0 K]);
% phi   = chebfun2(@(s,la) sum(a(k)*Cos(la*k)*j0(2.0*k*sqrt(s))), [0 Ls 0 La], 'vectorize');
% 
% 
%  figure(3);
%  plot(phi);
%  xlabel('$s$', 'interpreter', 'LaTeX', 'fontsize', 12);
%  ylabel('$\lambda$', 'interpreter', 'LaTeX', 'fontsize', 12);
%  view([0 90]); colorbar;
%  title('Phi Two-parameters integral $f(s,\lambda)$', 'interpreter', 'LaTeX', 'fontsize', 12);
%  set(gcf, 'color', 'w');
%  export_fig('phi.png', '-m2', '-a4', '-painters');
% 
%  disp('done')
% disp(phi(1,2))
% 
%  disp('psi...')
% 
% Sin = chebfun(@(lk) sin(lk), [0 La*K], 'vectorize');
% psi = chebfun2(@(s,la) s^(-1/2)*sum(a(k)*Sin(la*k)*j1(2.0*k*sqrt(s))), [0.01 Ls 0.01 La], 'vectorize');
% 
% disp('done')
% 
%  figure(4);
%  plot(psi);
%  xlabel('$s$', 'interpreter', 'LaTeX', 'fontsize', 12);
%  ylabel('$\lambda$', 'interpreter', 'LaTeX', 'fontsize', 12);
%  view([0 90]); colorbar;
%  title('Psi Two-parameters integral $f(s,\lambda)$', 'interpreter', 'LaTeX', 'fontsize', 12);
%  set(gcf, 'color', 'w');
%  export_fig('psi.png', '-m2', '-a4', '-painters');

save('psi_phi')


disp('data projection onto lamba = 0...')







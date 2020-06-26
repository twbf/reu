load('psi_phi_projection.mat')

len = 160;

g = sqrt(9.81);

x = zeros(len, len);
t = zeros(len, len);
h = zeros(len, len);
u = zeros(len, len);

lambda = linspace(0, 10, len);
s = linspace(0, 10, len);

u = psi(s,lamda);
h = phi(s,lamda)-u^2/2;

x = s - h;
t = (u + lamda)/g;

hh = reshape(h, [len*len,1]);
xx = reshape(x, [len*len,1]);
tt = reshape(t, [len*len,1]);
uu = reshape(u, [len*len,1]);




% for i = 1:len
%     for j = 1:len
%         
%         s = i/16 ;
%         
%         lamda = j/16 -1/16;
%         
%         u = psi(s,lamda);
% 
%         h = phi(s,lamda)-u^2/2;
% 
%         x = s - h;
% 
%         t = (u + lamda)/g;
%         
%         xx(i,j) = x;
%         tt(i,j) = t;
%         hh(i,j) = h;
%         uu(i,j) = u;
%         
%     end
% end
%scatter(xx,tt)

figure(1);
mesh(x,t,h);
xlabel('$k$', 'interpreter', 'LaTeX', 'fontsize', 12);
ylabel('$a(k)$', 'interpreter', 'LaTeX', 'fontsize', 12);
title('Coefficient $a(k)$', 'interpreter', 'LaTeX', 'fontsize', 13);
export_fig('Free Surface area', '-m2', '-a4', '-painters');

figure(2);
mesh(x,t,u);
xlabel('$k$', 'interpreter', 'LaTeX', 'fontsize', 12);
ylabel('$a(k)$', 'interpreter', 'LaTeX', 'fontsize', 12);
title('Coefficient $a(k)$', 'interpreter', 'LaTeX', 'fontsize', 13);
export_fig('Speed', '-m2', '-a4', '-painters');

figure(3);
scatter(xx,tt);
xlabel('$k$', 'interpreter', 'LaTeX', 'fontsize', 12);
ylabel('$a(k)$', 'interpreter', 'LaTeX', 'fontsize', 12);
title('Coefficient $a(k)$', 'interpreter', 'LaTeX', 'fontsize', 13);
export_fig('Free Surface area', '-m2', '-a4', '-painters');


ana = scatteredInterpolant(xx,tt,hh);
save('ana_interp1_0', 'ana')


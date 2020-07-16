%--------- Graphing Boundary Data Vector---------

close all
clear all
format Long

load('str.mat','shelf')

%----Paramters------%
g = 9.81;
g2 = 0.5*g;
d=1;

%---------Removing Dimensions------%
t_dim = linspace(0,10,200);
t = t_dim/(sqrt(d/g)); %dimensionless time

u1_dim = shelf.u1;                                  %pulls u1 from runwave.m
u1 = u1_dim/(sqrt(g*d)); %dimensionless velocity

eta1_dim = shelf.eta1; % pulls eta from runwave.m
eta1 = eta1_dim/(d); %making height dimensionless

lambda = t-u1; %lambda definition

%----------Some Graphs------------%

figure(1) %plot u
plot(t,u1)
grid on
title('Bounrdary speed data, u1(t)')
xlabel('t')
ylabel('u1(t)')

figure(2) %plot lambda
plot(t,lambda)
grid on
title('lambda(t) = t-u1(t)')
xlabel('t')
ylabel('lambda')  %lambda(t)

figure(3) %plot the inverse
plot(lambda,t)
grid on       %t(lambda), the inverse!
title('The inverse of lambda, t(lambda) (gamma)')
xlabel('lambda')
ylabel('t')

%----------Polynomial Fitting to find Gamma-------%

coeff = polyfit(lambda,t,1); % fitting a function to the inverse
c1 = coeff(1)
c2 = coeff(2)

gamma = c1*lambda+c2;  % gamma vector
gamma(1) = 0; %FIXME supposed to keep gama positive
gamma(200) = 31.163625990159890; %FIXME supposed to keep gama positive in range

%--------------Calculating Eta at Gamma-------%

eta_of_gamma = interp1(t,eta1,gamma); % calculating eta at gamma
disp(gamma);
figure(4)
plot(gamma,eta_of_gamma)
title('Eta of Gamma')
xlabel('Gamma')
ylabel('Eta')
grid on

sigma = eta_of_gamma+1;  %finding sigma

diff = gamma-lambda;
energy = eta_of_gamma+(0.5*diff.^2); % the "energy" term from the vector
psizero = [diff;energy];
squareterm = (0.5*lambda).^2;


%-----------Bonus: Graphing a few lines from the notes------%


% Trying to plot problematic area from "Ideas on Numerical" where Det(D)=0

% figure(5)         %What does this mean for us?? Not sure yet
% hold on
% plot(lambda,eta_of_gamma,'r')
% plot(lambda,lambda+squareterm,'b')
% title('Plotting Eta and Lambda+1/2Lambda^2')
% legend('Eta of Gamma','Lambda+1/2*lambda^2')
% xlim([-0.25,2.5])
% ylim([-10,15])
% xlabel('Lambda')
% hold off
% grid on

%---------Coding D--------%

n = length(t);
sigmaprime = zeros(1,200);
sigmaprime(1) = (sigma(2)-sigma(1))/(t(2)-t(1)); %initial point

for i = 2:n-1    %calclating sigma prime via finite difference

   sigmaprime(i) = (sigma(i+1)-sigma(i-1))/(t(i+1)-t(i-1));

end
sigmaprime(n) = (sigma(n)-sigma(n-1))/(t(n)-t(n-1)); %final point

figure(5)
plot(lambda,sigmaprime)
title('Eta of Gamma')
xlabel('Gamma')
ylabel('Eta')
grid on

minusone = repelem(-1,200); %vector full of -1

%----------Inverting D-------%

D_inverse = zeros(2,400); %preparing resultant D-inverse


for i = 1:200 % inverting element-by-element

D_sub = [sigmaprime(i),minusone(i);-(sigma(i)),sigmaprime(i)];

D_sub_inv = inv(D_sub);

D_inverse(1,i) = D_sub_inv(1,1);
D_inverse(1,200+i) = D_sub_inv(1,2); % because the resultant matrix is 2x400
D_inverse(2,i) = D_sub_inv(2,1);
D_inverse(2,200+i) = D_sub_inv(2,2);

end

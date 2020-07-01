

%initial conditions

H1 = 0.006;
H2 = 0.018;
c1 = 0.4444;
c2 = 4.0;
x1 = 4.1209;
x2 = 1.6384;

eta = chebfun(@(x) H1*exp(-c1*(x - x1).^2) - H2*exp(-c2*(x - x2).^2), [0 L]);
eta_prime = chebfun(@(x) -2*H1*c1*(x-x1)*exp(-c1*(x - x1).^2) + -2*H2*c2*(x-x2)*exp(-c2*(x - x2).^2), [0 L]);

%u   = chebfun(@(x) -eta(x)/sqrt(x), [0 L]);
%u_prime = chebfun(@(x) -( ( eta(x)/(2*sqrt(x)) + sqrt(x)*eta_prime(x) )/x ), [0 L]);  %needs to change

u = chebfun(@(x) 0, [0 L]);
u_prime = chebfun(@(x) 0, [0 L]);


%numeric solution

%Nicolsky 2018 soultion

%catalina 1 solution

%comparison

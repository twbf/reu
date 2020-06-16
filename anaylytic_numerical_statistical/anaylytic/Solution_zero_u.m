clear all

a=0.017;
b=4;
x_0=1.69;
m=1000;
c=sqrt(m/(m+1));                    %wave speed

atol = 1e-5;
rtol = 1e-5;
ntol = 1e-3;

slin=[0.001:.05:4];
taulin=[0.00:.05:0.5];



eta_0=@(x) a*exp(-b*(x-x_0).^2);
deta_0=@(x) -b*2*a.*(x-x_0).*exp(-b*(x-x_0).^2);

afun=@(x,k) 2*k*eta_0(x).*(x+eta_0(x)).^(1/(2*m)).*besselj(1/m,2*k*sqrt(x+eta_0(x))).*(1+deta_0(x));
%a = @(k) quadgk(@(x) afun(x,k), 0, 10, 'AbsTol', 1e-10, 'RelTol', 1e-10 );
a = @(kk) arrayfun( @(k) quadgk(@(x) afun(x,k), 0, 10, 'AbsTol', atol/10, 'RelTol', rtol/10), kk);

% %Plotting the coefficient "a(k)"
% kv=linspace(0,50,250);
% plot(kv, a(kv))


psi_fun=@(s,tau,k) s.^(-1/m/2).*besselj(1/m,2*k*sqrt(s)).*(a(k).*cos(c*k*tau));
psi = @(s,tau) quadgk(@(k) psi_fun(s,tau,k),  0, 50, 'AbsTol', atol, 'RelTol', rtol);

psi_s_fun=@(s,tau,k) -1/2*s.^(-(1/m+1)/2).*besselj(1/m+1,2*k*sqrt(s))*2.*k.*(a(k).*cos(c*k*tau));
psi_s = @(s,tau) quadgk(@(k) psi_s_fun(s,tau,k),  0, 50, 'AbsTol', atol, 'RelTol', rtol);


%phi_fun=@(s,tau,k) 1/c*s.^(-1/(2*m)-1/2).*besselj(1/m+1,2*k*sqrt(s)).*(a(k).*sin(c*k*tau));
phi_fun=@(s,tau,k) 1/c*s.^(-(1/m+1)/2).*besselj(1/m+1,2*k*sqrt(s)).*(a(k).*sin(c*k*tau));
phi = @(s,tau) quadgk(@(k) phi_fun(s,tau,k),  0, 50, 'AbsTol', atol, 'RelTol', rtol);

phi_tau_fun=@(s,tau,k) s.^(-(1/m+1)/2).*besselj(1/m+1,2*k*sqrt(s)).*(a(k).*cos(c*k*tau).*k);
phi_tau = @(s,tau) quadgk(@(k) phi_tau_fun(s,tau,k),  0, 50, 'AbsTol', atol, 'RelTol', rtol);

phi_s_fun=@(s,tau,k) -1/2/c*s.^(-(1/m+2)/2).*besselj(1/m+2,2*k*sqrt(s)).*2.*k.*(a(k).*sin(c*k*tau));
phi_s = @(s,tau) quadgk(@(k) phi_s_fun(s,tau,k),  0, 50, 'AbsTol', atol, 'RelTol', rtol);


% [sv,tauv]=meshgrid(slin, taulin);
% for i=1:size(sv,1)
%     for j=1:size(sv,2)
%         psiv(i,j)=psi(sv(i,j), tauv(i,j));
%         phiv(i,j)=phi(sv(i,j), tauv(i,j));
%     end
%     i
% end
% 
% uv=phiv;
% etav=psiv-uv.^2/2;
% xv=sv-etav;
% tv=tauv+uv;
% 
% 
% Eta = scatteredInterpolant(xv(:),tv(:),etav(:),'natural','none');
% U = scatteredInterpolant(xv(:),tv(:),uv(:),'natural','none');
% 
% xp = -0.1:.005:4;
% 
% trunup=0.47;
% eta_runup = Eta(xp,trunup*ones(size(xp)));
% u_runup = U(xp,trunup*ones(size(xp)));
% plot(xp, eta_runup)

tt=[1:0.025:6];

for ii=1:length(tt)
 
    tau=tt(ii); 
    phiv = phi(1e-10, tau);
    psiv = psi(1e-10, tau);

    tcoast(ii) = tau + phiv;
    xcoast(ii) = 0   - psiv + phiv^2/2;
    [tcoast(ii) xcoast(ii)]
end
% return
% save(['WaveData_coast.mat'],'xcoast','tcoast');



tt=[1.65 2.0 2.475 3.06];
for jj=1:length(tt)
    xx=linspace(spline(tcoast,xcoast,tt(jj)),10,1000);
    
    for ii=1:length(xx)
        xs=xx(ii);
        ts=tt(jj);
        
        if(ii==1), tau=ts; s=xs+1; t=100; x=100; end
        while ((abs(t-ts)>ntol)||(abs(x-xs)>ntol))
            phiv = phi(s, tau);
            psiv = psi(s, tau);
            
            t = tau + phiv;
            x = s   - psiv + phiv^2/2;
            
            tau = tau - (ts-t)/(-1-phi_tau(s,tau));
            s   = max(s   - (xs-x)/(-1+psi_s(s,tau)-phiv*phi_s(s,tau)), 1e-10);
            
            display(num2str([ii t ts x xs psiv-phiv^2/2 eta_0(x)],'%6.4f '))
        end
        phiv = phi(s, tau);
        psiv = psi(s, tau);
        
        uvv(ii) = phiv;
        etavv(ii) = psiv-phiv^2/2;
    end
    
    save(['WaveData_t=',num2str(tt(jj)),'.mat'],'xx','uvv','etavv');
end



clear tt xx uvv etavv;
xx=0.25;
tt=linspace(0,6,2000);

for ii=1:length(tt)
    xs=xx;
    ts=tt(ii);
    
    if(ii==1), tau=ts; s=xs+1; t=100; x=100; end
    while ((abs(t-ts)>ntol)||(abs(x-xs)>ntol))
        phiv = phi(s, tau);
        psiv = psi(s, tau);
        
        t = tau + phiv;
        x = s   - psiv + phiv^2/2;
        
        tau = tau - (ts-t)/(-1-phi_tau(s,tau));
        s   = max(s   - (xs-x)/(-1+psi_s(s,tau)-phiv*phi_s(s,tau)), 1e-10);
        
        display(num2str([ii t ts x xs psiv-phiv^2/2],'%6.4f '))
    end
    phiv = phi(s, tau);
    psiv = psi(s, tau);
    
    uvv(ii) = phiv;
    etavv(ii) = psiv-phiv^2/2;
end
save(['WaveData_x=',num2str(xx),'.mat'],'xx','tt','uvv','etavv');





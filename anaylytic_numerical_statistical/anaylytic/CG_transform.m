load('psi_phi_projection.mat')

len = 160;

g = sqrt(9.81);

xx = zeros(len, len);

tt = zeros(len, len);

hh = zeros(len, len);


for i = 1:len
    for j = 1:len
        
        s = i/16 ;
        
        lamda = j/16 -1/16;
        
        u = psi(s,lamda);

        h = phi(s,lamda)-u^2/2;
        
        if h <= -1
            h = -1;
        end

        x = s - h;

        t = (u + lamda)/g;
        
        xx(i,j) = x;
        tt(i,j) = t;
        hh(i,j) = h;
        uu(i,j) = u;
        
    end
end
%scatter(xx,tt)



mesh(xx,tt,hh)

title(['$\eta(x,t)$ by CG Transform'], IN, 'latex', FS, 14);
xlabel('$x$', IN, 'latex', 'fontsize', 16);
ylabel('$t$', IN, 'latex', 'fontsize', 16);


hh = reshape(hh, [len*len,1]);

xx = reshape(xx, [len*len,1]);

tt = reshape(tt, [len*len,1]);

ana = scatteredInterpolant(xx,tt,hh,'natural','none')
save('ana_interp1_0', 'ana')

view(3)


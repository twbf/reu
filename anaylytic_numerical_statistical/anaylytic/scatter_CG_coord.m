load('psi_phi_old.mat')

len = 160;

xx = zeros(len, len);

tt = zeros(len, len);

hh = zeros(len, len);


for i = 1:len
    for j = 1:len
        
        s = i/16;
        
        lamda = j/16;
        
        u = psi(s,lamda);

        h = phi(s,lamda)-u^2/2;
        
        if h <= -1
            h = -1;
        end

        x = s - h;

        t = u + lamda;
        
        xx(i,j) = x;
        tt(i,j) = t;
        hh(i,j) = h;
        
    end
end
% 
hh = reshape(hh, [len*len,1]);

xx = reshape(xx, [len*len,1]);

tt = reshape(tt, [len*len,1]);

ana = scatteredInterpolant(xx,tt,hh)
save('ana_interp', 'ana')

scatter(xx,tt)


title(['CG Transformed Coordinates $10 X 10$ Grid in $(s,\lambda)$ with $1/4$ step sizeto $(x,t)$'], IN, 'latex', FS, 14);
xlabel('$x$', IN, 'latex', 'fontsize', 16);
ylabel('$t$', IN, 'latex', 'fontsize', 16);

%x = gallery('uniformdata',[20,1],0);
%y = gallery('uniformdata',[20,1],1);
DT = delaunay(xx,tt);
trimesh(DT,xx,tt,hh);


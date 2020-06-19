load('ana_interp1')
load('num_interp1')

size = 1000;

dx = 1/110;
dt = 1/330;

diff = zeros(size,size);

l2 = zeros(size);

numerical = zeros(size,size);

anaylytic = zeros(size,size);

for j=1:size
    for i=1:size
    
        x = dx*i;
        t = dt*j;
        
        numerical(i,j) = num(x-2,t);
        
        if numerical(i,j) == 0
            anaylytic(i,j) == 0;
        else
            anaylytic(i,j) = ana(x-2.5,t);
        end
        
        diff(i,j) = anaylytic(i,j)-numerical(i,j);
        
    end
    
    l2(j) = norm(diff(:,j),2);
end


xx = zeros(size,size);

tt = zeros(size,size);

tt_1 = zeros(size);

for i=1:size
    tt(i,:) = linspace(0,size*dt,size);
    
    tt_1 = linspace(0,size*dt,size);
end

for j=1:size
    xx(:,j)= linspace(-2.5,size*dx-2.5,size);
end




figure(1)
mesh(xx,tt,numerical)
title(['Numerical Solution (Deny FV)'], IN, 'latex', FS, 14);
xlabel('$x$', IN, 'latex', 'fontsize', 16);
ylabel('$t$', IN, 'latex', 'fontsize', 16);

figure(2)
mesh(xx,tt,anaylytic)
title(['Anaylytical Solution (Nicolsky et al. 2018)'], IN, 'latex', FS, 14);
xlabel('$x$', IN, 'latex', 'fontsize', 16);
ylabel('$t$', IN, 'latex', 'fontsize', 16);

figure(3)
mesh(xx,tt,diff)
title(['Anaylytical-Numerical Difference'], IN, 'latex', FS, 14);
xlabel('$x$', IN, 'latex', 'fontsize', 16);
ylabel('$t$', IN, 'latex', 'fontsize', 16);


figure(4)
plot(tt_1, l2)
title(['L2 Norm of the differnce for each $t$'], IN, 'latex', FS, 14);
xlabel('$t$', IN, 'latex', 'fontsize', 16);
ylabel('L2 Norm', IN, 'latex', 'fontsize', 16);






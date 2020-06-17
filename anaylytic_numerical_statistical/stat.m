load('ana_interp')
load('num_interp')

size = 100;

dx = 1/10;

diff = zeros(size,size);

for i=1:size
    for j=1:size
        x = dx*i;
        t = dx*j;
        
        diff(i,j) = ana(x,t)-num(x,t);
        
    end
end
        
mesh(diff)

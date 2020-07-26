t0 = 0.0;
Tf = 10.0;
M     = 2000;	% number of time instances where we project solution
tlist = linspace(t0, Tf, M);

shelf.u1 = sin(tlist);
shelf.eta1 = 0.004*cos(tlist)+1;

save('str_sine.mat','shelf')

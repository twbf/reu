%{
    Finite volume solver for NSWE wave propagation and run-up
    Copyright (C) 2020 Denys DUTYKH

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

%}

%%% -------------------------------------------------- %%%
%%% The main file for NSWE solver based on FV method   %%%
%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, CNRS -- LAMA, Univ of Savoie %%%
%%% E-mail: Denys.Dutykh@univ-savoie.fr                %%%
%%% Web:    http://www.denys-dutykh.com/               %%%
%%% Blog:   http://dutykh.github.io/                   %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%
%%% Distributed under GNU General Public License       %%%
%%% -------------------------------------------------- %%%

clear
close all
format longE

%%% Libraries we use:
addpath('sources/');
addpath('odetpbar/');
addpath('export_fig/');

%%% Global variables:
global cf2 d xc FS IN LW
global a amp td b g g2 dx h N

FS = 'FontSize';
IN = 'Interpreter';
LS = 'LineStyle';
LW = 'LineWidth';

%%% Physical parameters:
g  = 1.0;	% gravity acceleration
g2 = 0.5*g;	% g/2
cf2 = 0.0;	% friction coefficient

a  = -1.0;	% the left boundary (incident wave)
b  = 0.45;	% the right boundary (wall)
d  = 1.0;   % still water depth on the right boundary
% Bottom slope:
td = 1.0/1.0;

%%% Numerical parameters:
N  = 1000;							% number of grid points
x  = linspace(a, b, N+1)';			% cell interfaces
dx = x(2) - x(1);                  	% spatial grid step
xc = 0.5*(x(1:end-1) + x(2:end));  	% centers of cells

%%% Bathymetry function:
h  = -td*xc;

%%% Choice of the initial condition:
w0 = zeros(2*N,1);
w0(1:N) = max(h, eps+0*h);   % zero initial condition without velocities

% time stepping:
t0 = 0.0;
Tf = 10.0; % final simulation time (the same as experiment)

%%% Plot the initial condition to check:
amp = 0.1;
figure;
set(gcf, 'pos', [1 621 903 353]);
Plot(t0, w0);

%%% We run the simulation:
options = odeset('AbsTol', 1e-4, 'RelTol', 1e-4, 'OutputFcn', @odetpbar, 'MaxStep', 1.0);
fprintf('Simulation...\n');
sol = ode23(@RHS, [t0 Tf], w0, options);
fprintf(' Done\n');

%%% Post-processing of the solution:
M     = 500;	% number of time instances where we project solution
tlist = linspace(t0, Tf, M);
solpr = deval(sol, tlist);

% %%% Make a little animation of the numerical solution:
Rup = zeros(M,1);
for t=1:M % loop in time
    ind = find(solpr(1:N,t) > 1e-2, 1, 'last');
    Rup(t) = -h(ind);
	Plot(tlist(t), solpr(:,t));
end % for t
Rup = Rup - Rup(1);

%%% Extraction of run-up data:
Rup = smooth(Rup, 5);
figure;
plot(tlist, Rup, '-', LW, 2.0); grid off;
axis tight;
xlim([t0 Tf]); ylim([1.2*min(Rup) 1.1*max(Rup)]);
xlabel('$t$, s', IN, 'LaTeX', FS, 14);
ylabel('$R(t)$, m', IN, 'LaTeX', FS, 14);
ll = legend(' Numerical data', 'location', 'NorthWest');
set(ll, 'box', 'off');
title('Wave run-up');
set(gcf, 'Color', 'w');

% export_fig('shots/runup.png', '-m2', '-a4', '-painters', '-q100');
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

function [eta, u] = run_num()

    %%% Libraries we use:
    addpath('IV_FVM/sources/');
    addpath('IV_FVM/odetpbar/');

    %%% Global variables:
    %for numerical scheme
    global cf2 d xc FS IN LW
    global a amp td b g g2 dx h N

    %IC variables from run.m
    global eta_0 eta_prime u_0 u_prime td t0 Tf x0 Xf
    global x_res t_res %resolution

    fprintf('Numeric Simulation...\n');

    FS = 'FontSize';
    IN = 'Interpreter';
    LS = 'LineStyle';
    LW = 'LineWidth';

    %%% Physical parameters:
    g2 = 0.5*g;	% g/2
    cf2 = 0.0;	% friction coefficient

    %a  = -1;	% the left boundary (incident wave)
    %b  = 10;	% the right boundary (wall)

    %%% Numerical parameters:
    N  = x_res;							% number of grid points
    x  = linspace(x0, Xf, N+1)';			% cell interfaces (the apostrophe is to transpose)
    dx = x(2) - x(1);                  	% spatial grid step
    xc = 0.5*(x(1:end-1) + x(2:end));  	% centers of cells

    %%% Bathymetry function:
    h  = td*xc;

    mm=x_res;

    small_h = h(1:mm);

    %%% Choice of the initial condition:
    w0 = zeros(2*N,1);

    w0(1:N) = max(h, eps+0*h);   % zero initial condition without velocities

    w0(1:N) = w0(1:N) + eta_0(xc);

    %setting speed
    %u = 0 where x<0
    for i = 1:N
      if xc(i) > 0
        w0(i+N) = w0(i)*u_0(xc(i));
      end
    end

    %w0(N+1:2*N,1) = w0(1:N).*u_0(xc);

    %%% Plot the initial condition to check:
    %amp = 1.5*H2;
    %figure;
    %set(gcf, 'pos', [1 621 903 353]);
    %Plot(t0, w0);

    %%% We run the simulation:
    options = odeset('AbsTol', 1e-4, 'RelTol', 1e-4, 'OutputFcn', @odetpbar, 'MaxStep', 1.0);
    sol = ode23(@RHS, [t0 Tf], w0, options);

    %%% Post-processing of the solution:
    M     = t_res;	% number of time instances where we project solution
    tlist = linspace(t0, Tf, M);
    solpr = deval(sol, tlist);

    eta = zeros(mm,M);
    u = zeros(mm,M);
    for i=1:mm
      for j=1:M

        u(i,j) = solpr(N + i,j);

        if solpr(i,j)<0.0000001 &&  solpr(i,j)>-0.0000001
            eta(i,j) = 0;
        else
            eta(i,j) = (solpr(i,j)-small_h(i));
        end
      end
    end

end

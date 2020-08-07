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
    addpath('sources/');
    addpath('odetpbar/');

    %%% Global variables:
    %for numerical scheme
    global cf2 d xc FS IN LW
    global a amp td b g g2 dx h N
    %IC variables from run.m
    global H1 H2 c1 c2 x1 x2 eta_0 u_0 td t0 Tf x0 Xf

    eta_0 = @(x) 0.1*exp(-(x-5).^2);
    u_0 = @(x) 0;
    td = 1;
    t0 = 0;
    Tf = 3.1;



    H1 = 0.006;
    H2 = 0.018;
    c1 = 0.4444;
    c2 = 4.0;
    x1 = 4.1209;
    x2 = 1.6384;

    eta_0 = chebfun(@(x) H1*exp(-c1*(x - x1).^2) - H2*exp(-c2*(x - x2).^2), [-2, 20]);
    eta_prime = diff(eta_0);

    disp(" ")
    disp("Deny's Numeric:");

    FS = 'FontSize';
    IN = 'Interpreter';
    LS = 'LineStyle';
    LW = 'LineWidth';

    %%% Physical parameters:
    g  = 9.8;	% gravity acceleration
    g2 = 0.5*g;	% g/2
    cf2 = 0.0;	% friction coefficient

    a  = -1;	% the left boundary (incident wave)
    b  = 10;	% the right boundary (wall)

    %%% Numerical parameters:
    N  = 8000;							% number of grid points
    x  = linspace(a, b, N+1)';			% cell interfaces (the apostrophe is to transpose)
    dx = x(2) - x(1);                  	% spatial grid step
    xc = 0.5*(x(1:end-1) + x(2:end));  	% centers of cells

    %%% Bathymetry function:
    h  = td*xc;

    mm=8000;

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
    amp = 1.5*H2;
    %figure;
    %set(gcf, 'pos', [1 621 903 353]);
    %Plot(t0, w0);

    %%% We run the simulation:
    options = odeset('AbsTol', 1e-4, 'RelTol', 1e-4, 'OutputFcn', @odetpbar, 'MaxStep', 1.0);
    fprintf('Numeric Simulation...\n');
    sol = ode23(@RHS, [t0 Tf], w0, options);
    fprintf(' Done\n');

    %%% Post-processing of the solution:
    M     = 800;	% number of time instances where we project solution
    tlist = linspace(t0, Tf, M);
    solpr = deval(sol, tlist);

    % %%% Make a little animation of the numerical solution:
    %Rup = zeros(M,1);
    %for t=1:M % loop in time
      %ind = find(solpr(1:N,t) > 1e-2, 1, 'first');

      %Rup(t) = -h(ind);

	     %Plot(tlist(t), solpr(:,t));
    %end % for t

    %Rup = Rup - Rup(1);


    x_mat = zeros(mm,M);

    t_mat = zeros(mm,M);

    for i=1:mm
      t_mat(i,:) = tlist;
    end

    for j=1:M
      x_mat(:,j)= linspace(-2, 10, mm);
    end

    %adding bathymetry

    %figure(1);

    hh = zeros(mm,M);
    uu = zeros(mm,M);
    for i=1:mm
      for j=1:M

        uu(i,j) = solpr(N + i,j);

        if solpr(i,j)<0.0001 &&  solpr(i,j)>-0.0001
            hh(i,j) = 0;
        else
            hh(i,j) = (solpr(i,j)-small_h(i));
        end
      end
    end

    eta_fvm = hh;
    save('fvm_test', 'eta_fvm');

    %     figure(2);
    %     mesh(hh);
    %     title(['$\eta(x,t)$ by FVM'], IN, 'latex', FS, 14);
    %     xlabel('$x$', IN, 'latex', 'fontsize', 16);
    %     ylabel('$t$', IN, 'latex', 'fontsize', 16);

    figure(1);
    mesh(hh);

    hh = reshape(hh, [mm*M,1]);
    uu = reshape(uu, [mm*M,1]);
    xx = reshape(x_mat, [mm*M,1]);
    tt = reshape(t_mat, [mm*M,1]);

    eta = scatteredInterpolant(xx,tt,hh);

    u = scatteredInterpolant(xx,tt,uu);
    % save('num_interp_cat1_0u_1s', 'num')

    %%% Extraction of run-up data:
    % Rup = smooth(Rup, 7);
    % figure;
    % plot(tlist, Rup, '-', LW, 2.0); grid off;
    % axis tight;
    % xlim([t0 Tf]); ylim([1.2*min(Rup) 1.1*max(Rup)]);
    % xlabel('$t$, s', IN, 'LaTeX', FS, 14);
    % ylabel('$R(t)$, m', IN, 'LaTeX', FS, 14);
    % ll = legend(' Numerical data', 'location', 'SouthEast');
    % set(ll, 'box', 'off');
    % title('Wave run-up');
    % set(gcf, 'Color', 'w');

end

%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, CNRS -- LAMA, Univ of Savoie %%%
%%% E-mail: Denys.Dutykh@univ-savoie.fr                %%%
%%% Web:    http://www.denys-dutykh.com/               %%%
%%% Blog:   http://dutykh.github.io/                   %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%
%%% Distributed under GNU General Public License       %%%
%%% -------------------------------------------------- %%%

function Plot (t, w)

    global FS IN LW	
	global a amp b td h xc N

	plot(xc, w(N:2*N-1) - h, '-', LW, 2.2), hold on, grid on
    plot(xc, -h, '-', LW, 2.2), hold off
    xlim([-0.5 7]);
    xlabel('$x$', IN, 'latex', 'fontsize', 16);
    ylabel('$\eta(x,t)$', IN, 'latex', FS, 16);
    title(['Free surface elevation at $t = $ ',num2str(t,'%5.2f'), ' s'], IN, 'latex', FS, 14);
    set(gcf, 'Color', 'w');
    drawnow

end % Plot ()
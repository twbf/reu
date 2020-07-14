function Plot_xt(f, size, figure_num, graph_title, x_axis, y_axis, z_axis)

    global t0 Tf x0 Xf

    num_x = size;
    num_t = size;

    t_list = linspace(t0, Tf, num_t);
    x_list = linspace(x0, Xf, num_x);

    xx = zeros(num_x, num_t);
    tt = zeros(num_x, num_t);
    hh = zeros(num_x, num_t);

    for t=1:num_t
      for x=1:num_x
        xx(x,t) = t_list(t);
        tt(x,t) = x_list(x);
        hh(x,t) = f(x_list(x),t_list(t));
      end
    end

    figure(figure_num)
    mesh(tt,xx,hh)
    title([graph_title],'Interpreter','latex');
    xlabel(x_axis,'Interpreter','latex');
    ylabel(y_axis,'Interpreter','latex');
    zlabel(z_axis,'Interpreter','latex');

end

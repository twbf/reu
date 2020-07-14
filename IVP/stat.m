
function [diff, l2, g_diff] = stat(f, g, size, nuclear_option)

    global t0 Tf x0 Xf

    num_x = size;
    num_t = size;

    t_list = linspace(t0, Tf, num_t);
    x_list = linspace(x0, Xf, num_x);

    difference = zeros(size,size);
    l2 = zeros(size);
    g_adj = zeros(size,size);

    for t=1:num_t
        for x=1:num_x
          if nuclear_option
            if f(x_list(x),t_list(t))  == 0
              difference(x,t) = 0;
              g_adj(x,t) = 0;
            else
              g_adj(x,t) = g(x_list(x),t_list(t));
              difference(x,t) = f(x_list(x),t_list(t)) - g(x_list(x),t_list(t));
            end
          else
            difference(x,t) = f(x_list(x),t_list(t)) - g(x_list(x),t_list(t));
          end
        end
        l2(t) = norm(difference(:,t),2);
    end

    xx = zeros(size, size);
    tt = zeros(size, size);

    for i=1:size
      tt(i,:) = t_list;
      xx(:,i) = x_list;
    end

    xx = reshape(xx, [size * size,1]);
    tt = reshape(tt, [size * size,1]);
    difference = reshape(difference, [size * size,1]);
    g_adj = reshape(g_adj, [size * size,1]);

    diff = scatteredInterpolant(xx,tt,difference);
    g_diff = scatteredInterpolant(xx,tt,g_adj);

end

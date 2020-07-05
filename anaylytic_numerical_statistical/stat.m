
function [diff, l2] = stat(f, g, size)

    global t0 Tf x0 Xf

    num_x = size;
    num_t = size;

    t_list = linspace(t0, Tf, num_t);
    x_list = linspace(x0, Xf, num_x);

    difference = zeros(size,size);
    l2 = zeros(size);

    for t=1:num_t
        for x=1:num_x
            difference(x,t) = f(x_list(x),t_list(t)) - g(x_list(x),t_list(t));
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

    diff = scatteredInterpolant(xx,tt,difference);

end

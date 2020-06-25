% proj[1] = psi    proj[2] = phi

function phi = phi_0 (x)

    % input array x is a list of values for the projection 
    % note: memory should be reaalocated

    global proj

    for i=1:size(x,2)
        tmp = proj(x(i));
        phi(i) = tmp(1);
    end
    
end
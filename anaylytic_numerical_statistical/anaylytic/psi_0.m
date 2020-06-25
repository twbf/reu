function psi = psi_0 (x)

    % input array x is a list of values for the projection 
    % note: memory should be reaalocated

    global proj

    for i=1:size(x,2)
        tmp = proj(x(i));
        psi(i) = tmp(2);
    end
    
end
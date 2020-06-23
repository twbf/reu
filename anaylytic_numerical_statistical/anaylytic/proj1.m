% proj[1] = psi    proj[2] = phi

function proj = proj1 (xx)

    global eta eta_prime u u_prime beta

    %s = xx + eta(xx);
    s = chebfun(@(x) x + eta(x), [0 20]);
    A = @(x) [0 1; beta^2*s(x) 0]; %beta^2*s
    B = [0 0; 1 0];
    
    uiu = chebfun(@(x) 1, [0 20]);
    
    disp("D")

    
    
    D = chebfun(@(x) ([1 0].*eta_prime(x))*eye(2) + ([1 0].*u_prime(x))*A(x), [0 20]);
    %D = chebfun(@(x) ([1 0].*eta_prime(x))*eye(2), [0 20]);
    %D = (1+eta_prime(2))*eye(2);
    disp(D)
    disp(D(:,:))
    
    disp(inv(D))
    phi0 = [u(xx) ; eta(xx)+(u(xx)^2)/2];
    phi0_prime = [u_prime(xx); eta_prime(xx)+2*u(xx)*u_prime(xx)];
    
    proj = phi0 + u(xx)*(u_prime(xx)*inv(D)*B*phi0 - B*phi0-A*inv(D)*phi0_prime);
    
end
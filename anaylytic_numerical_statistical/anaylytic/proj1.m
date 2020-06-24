% proj[1] = psi    proj[2] = phi

function proj = proj1 (xx)

    global eta eta_prime u u_prime beta

    s = chebfun(@(x) x + eta(x), [0 20]);
    A = @(x) [0 1; beta^2*s(x) 0]; %beta^2*s
    B = [0 0; 1 0];
    
    %D = chebfun(@(x) ([1 0].*eta_prime(x))*eye(2) + ([1 0].*u_prime(x))*A(x), [0 20]);
    
    D = @(x) eta_prime(x)*eye(2) + u_prime(x)*A(x);
    
    %D = chebfun(@(x) ([1 0].*eta_prime(x))*eye(2), [0 20]);
    %D = (1+eta_prime(2))*eye(2);
    %disp(D(1))
    %disp(inv(D(1)))
    
    phi0 = @(x) [u(x) ; eta(x)+(u(x)^2)/2];
    phi0_prime = @(x) [u_prime(x); eta_prime(x)+2*u(x)*u_prime(x)];
    
    proj = @(x) phi0(x) + u(x)*(u_prime(x)*inv(D(x)).*B.*phi0(x) - B.*phi0(x)-A.*inv(D(x)).*phi0_prime(x));
    
end
function [u] = U0(x)
    nu = 0.5;
    sigma = 0.1;
    u = 0.5/(sigma*sqrt(2*pi))*exp(-0.5*((x-nu)/sigma)^2);
    u = u+1;
end
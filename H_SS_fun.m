function[H_SS] = H_SS_fun(p,n)
    %p is the input parameters, and s is the section of the blade.
    B_SS = B_SS_fun(p,n);
    Theta_minus = Theta_minus_fun(p,n);
    k_bar_x = k_bar_x_fun(p,n);
    kappa_bar = kappa_bar_fun(p,n);
    alpha_SS = p.V./p.U_c;
    
    H_SS = @(freq) ((1+1i).*exp(-4.*1i.*kappa_bar).*(1-Theta_minus.^2))./(2.*sqrt(pi).*...
    (alpha_SS-1).*p.k_bar_x.*sqrt(B_SS));%k_x might be a vector; %Maybe change M
end
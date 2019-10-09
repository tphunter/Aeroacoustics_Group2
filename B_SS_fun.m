function[B_SS] = B_SS_fun(p,n)
    %p is the input parameters, and s is the section of the blade.
    K_x = K_x_fun(p,n);
    mu_bar = mu_bar_fun(p,n);
    kappa_bar = kappa_bar_fun(p,n);
    B_SS = @(freq) K_x - p.Mach*mu_bar + kappa_bar; %Maybe change
end
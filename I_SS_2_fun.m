function[I_SS_2] = I_SS_2_fun(p,n)
    %p is the input parameters, and s is the section of the blade.
    H_SS = H_SS_fun(p,n);
    G_SS = G_SS_fun(p,n);
    Theta_minus = Theta_minus_fun(p,n);
    k_bar_x = k_bar_x_fun(p,n);
    k_bar = k_bar_fun(p,n);
    mu_bar = mu_bar_fun(p,n);
    Corr_epsilon = Corr_epsilon_fun(p,n);
    
    I_SS_2 = @(freq) H_SS.*Corr_epsilon+H_SS.*(exp(-2.*1i.*theta)...
        +1i.*(Theta_minus+k_bar_x+p.Mach.*mu_bar-k_bar).*G_SS); 
    %Take a look at definition of epsilon (slide 8 lecture 4), it should not be an exponent.%k_x might be a vector; %Maybe change M
end
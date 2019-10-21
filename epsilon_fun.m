function[epsilon] = epsilon_fun(p,n,freq,psi)
    %p is the input parameters, and s is the section of the blade.
    kappa_bar = kappa_bar_fun(p,n,freq,psi);
    
    epsilon =1/sqrt(1+1/(4*kappa_bar));%k_x might be a vector; %Maybe change M
end
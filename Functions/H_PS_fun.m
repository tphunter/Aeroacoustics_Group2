function[H_PS] = H_PS_fun(p,n,freq,psi)
    %p is the input parameters, and s is the section of the blade.
    B_PS = B_PS_fun(p,n,freq,psi);
    Theta_minus = Theta_minus_fun(p,n,freq,psi);
    k_bar_x = k_bar_x_fun(p,n,freq,psi);
    kappa_bar = kappa_bar_fun(p,n,freq,psi);
    alpha_PS = p.V./p.U_c;
    
    H_PS = ((1+1i).*exp(-4.*1i.*kappa_bar).*(1-Theta_minus.^2))./(2.*sqrt(pi).*...
    (alpha_PS-1).*k_bar_x.*sqrt(B_PS));%k_x might be a vector; %Maybe change M
end
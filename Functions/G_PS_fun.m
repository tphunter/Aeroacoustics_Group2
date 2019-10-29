function[G_PS] = G_PS_fun(p,n,freq,psi)
    %p is the input parameters, and s is the section of the blade.
    Theta_minus = Theta_minus_fun(p,n,freq,psi);
    k_bar = k_bar_fun(p,n,freq,psi);
    kappa_bar = kappa_bar_fun(p,n,freq,psi);
    epsilon=epsilon_fun(p,n,freq,psi);
    
    G_PS = (1+epsilon).*exp(1i.*(2.*kappa_bar+Theta_minus)).*sin(Theta_minus-2.*kappa_bar)./...
    (Theta_minus-2.*kappa_bar)+(1-epsilon).*exp(1i.*(-2.*kappa_bar+Theta_minus)).*...
    sin(Theta_minus+2.*kappa_bar)./(Theta_minus+2.*kappa_bar)+(1+epsilon).*(1-1i)./...
    (2.*(Theta_minus-2.*kappa_bar)).*exp(4.*1i.*k_bar).*fresnel(4.*kappa_bar)-(1-epsilon)...
    .*(1+1i)./(2.*(Theta_minus+2.*kappa_bar)).*exp(4.*1i.*k_bar).*fresnel(4.*kappa_bar)...
    +exp(2.*1i.*Theta_minus)./2.*sqrt(2.*kappa_bar./Theta_minus).*fresnel(2.*Theta_minus).*...
    ((1-epsilon).*(1+1i)./(Theta_minus+2.*kappa_bar)-(1+epsilon).*(1-1i)./(Theta_minus-2.*kappa_bar));%k_x might be a vector; %Maybe change M
end
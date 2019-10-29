function[Theta_minus] = Theta_minus_fun(p,n,freq,psi)
    %p is the input parameters, and s is the section of the blade.
    kappa_bar = kappa_bar_fun(p,n,freq,psi);
    mu_bar = mu_bar_fun(p,n,freq,psi);
    sigma = sqrt(p.x.^2+(1-p.Mach.^2).*(p.y.^2+p.z.^2)); %Maybe change M
    Theta_minus = kappa_bar - mu_bar.*(p.x./sigma); %Maybe change M
end
function[C_SS] = C_SS_fun(p,n,freq,psi)
    %p is the input parameters, and s is the section of the blade.
    K_x = K_x_fun(p,n,freq,psi);
    mu_bar = mu_bar_fun(p,n,freq,psi);
    sigma=sqrt(p.x.^2+(1-p.Mach.^2).*(p.y.^2+p.z.^2)); %Maybe change M
    C_SS = K_x - mu_bar.*(p.x./sigma-p.Mach); %Maybe change M
end
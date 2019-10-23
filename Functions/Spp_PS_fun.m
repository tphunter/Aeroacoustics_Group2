function[Spp_PS] = Spp_PS_fun(p,n,freq,psi)
    
    sigma=sqrt(p.x.^2+(1-p.Mach.^2).*(p.y.^2+p.z.^2)); %Mach
    Phi_tot_PP_PS = Phi_tot_PP_PS_fun(p,n,freq,psi);
    k_bar_y = k_bar_y_fun(p,n,freq,psi);
    k_bar = k_bar_fun(p,n,freq,psi);
    I_PS_tot = I_PS_tot_fun(p,n,freq,psi);
    
    intS_PP_TE_PS= @(freq) (Phi_tot_PP_PS(n).*sin(p.R1./p.chord(n).*(k_bar_y-k_bar.*p.z./sigma)).^2./...
    ((k_bar_y-k_bar.*p.z./sigma).^2).*abs(I_PS_tot).^2);

    Spp_PS(n) = (k_bar(n).*p.z./(2.*pi.*sigma.^2)).^2.*2.*p.chord(n).*...
    integral(intS_PP_TE_PS,-100,100,'ArrayValued',1);
    
end
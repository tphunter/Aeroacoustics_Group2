function[Spp_SS] = Spp_SS_fun(p,n,freq,psi)
    
    k_bar(n) = k_bar_fun(p,n,freq,psi);
    sigma=sqrt(p.x.^2+(1-p.Mach.^2).*(p.y.^2+p.z.^2)); %Mach
    Phi_tot_PP_SS = Phi_tot_PP_SS_fun(p,n,freq,psi);
    k_bar_y = k_bar_y_fun(p,n,freq,psi);
    I_SS_tot = I_SS_tot_fun(p,n,freq,psi);
    
    intS_PP_TE_SS= @(freq) (Phi_tot_PP_SS(n).*sin(p.R1./p.chord(n).*(k_bar_y-k_bar.*p.z./sigma)).^2./...
    ((k_bar_y-k_bar.*p.z./sigma).^2).*abs(I_SS_tot).^2);

    Spp_SS = (k_bar(n).*p.z./(2.*pi.*sigma.^2)).^2.*2.*p.chord(n).*...
    integral(intS_PP_TE_SS,-1000,1000,'ArrayValued',1);
    
end
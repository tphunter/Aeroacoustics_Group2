function[Spp_PS] = Spp_PS_fun(p,n)
    
    k_bar = k_bar_fun(p,n);
    sigma=sqrt(p.x.^2+(1-p.Mach.^2).*(p.y.^2+p.z.^2)); %Mach
    Phi_tot_PP_PS = Phi_tot_PP_PS_fun(p,n);
    k_bar_y = k_bar_y_fun(p,n);
    k_bar = k_bar_fun(p,n);
    I_PS_tot = I_PS_tot_fun(p,n);
    
    intS_PP_TE_PS= @(freq) (Phi_tot_PP_PS.*sin(p.R1./p.chord(n).*(k_bar_y-k_bar.*p.z./sigma)).^2./...
    ((k_bar_y-k_bar.*p.z./sigma).^2).*abs(I_PS_tot).^2);

    Spp_PS = @(freq) (k_bar.*p.z./(2.*pi.*sigma.^2)).^2.*2.*p.chord(n).*...
    integral(intS_PP_TE_PS,-inf,inf,'ArrayValued',1); 
    
end
function[Phi_tot_PP_PS] = Phi_tot_PP_PS_fun(p,n,freq,psi)
    
    phi_PP_PS = phi_PP_PS_fun(p,n,freq,psi);
    l_y = l_y_fun(p,n,freq,psi);

    Phi_tot_PP_PS = phi_PP_PS.*l_y./pi;
end
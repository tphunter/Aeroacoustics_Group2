function[Phi_tot_PP_SS] = Phi_tot_PP_SS_fun(p,n,freq,psi)
    
    phi_PP_SS = phi_PP_SS_fun(p,n,freq,psi);
    l_y = l_y_fun(p,n,freq,psi);

    Phi_tot_PP_SS = phi_PP_SS.*l_y./pi;
end
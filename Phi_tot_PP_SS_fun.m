function[Phi_tot_PP_SS] = Phi_tot_PP_SS_fun(p,n)
    
    phi_PP_SS = phi_PP_SS_fun(p,n);
    l_y = l_y_fun(p,n);

    Phi_tot_PP_SS = @(freq) phi_PP_SS.*l_y./pi;
end
function[I_SS_tot] = I_SS_tot_fun(p,n,freq,psi)
    %p is the input parameters, and s is the section of the blade.
    I_SS_1 = I_SS_1_fun(p,n,freq,psi);
    I_SS_2 = I_SS_2_fun(p,n,freq,psi);
    
    I_SS_tot = I_SS_1 + I_SS_2;
end
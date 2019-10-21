function[I_PS_tot] = I_PS_tot_fun(p,n,freq,psi)
    %p is the input parameters, and s is the section of the blade.
    I_PS_1 = I_PS_1_fun(p,n,freq,psi);
    I_PS_2 = I_PS_2_fun(p,n,freq,psi);
    
    I_PS_tot = I_PS_1 + I_PS_2;
end
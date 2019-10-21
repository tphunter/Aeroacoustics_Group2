function[K_x] = K_x_fun(p,n,freq,psi)
    %p is the input parameters, and s is the section of the blade.
    K_x = freq*p.chord(n)/(2*p.U_c(n));
end
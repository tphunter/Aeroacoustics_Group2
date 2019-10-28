function[u_tau_PS] = u_tau_PS(p,n,freq,psi)
    %p is the input parameters, and s is the section of the blade.
    
    u_tau_PS = sqrt(p.tau_w_SS(n)/p.density);
end
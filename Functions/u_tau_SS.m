function[u_tau_SS] = u_tau_SS(p,n,freq,psi)
    %p is the input parameters, and s is the section of the blade.
    
    u_tau_SS = sqrt(p.tau_w_SS(n)/p.density);
end
function[k_bar_x] = k_bar_x_fun(p,n,freq,psi)
    k_bar_x = sin(p.phi)*cos(p.theta)*freq*p.c_R(n)*p.R1/(2*p.c);
end
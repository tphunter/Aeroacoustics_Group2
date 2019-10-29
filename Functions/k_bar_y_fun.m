function[k_bar_y] = k_bar_y_fun(p,n,freq,psi)
    k_bar_y = sin(p.phi)*sin(p.theta)*freq*p.c_R(n)*p.R1/(2*p.c);
end
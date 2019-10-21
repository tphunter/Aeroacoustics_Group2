function[PI_SS] = PI_SS_fun(p,n,freq,psi)

Karman = 0.41;

PI_SS=@(PI_SS_var) 2.*PI_SS_var-log(1+PI_SS_var)-Karman.*p.V./...
    p.u_tau-log(p.deltastar_SS.*p.V./p.kinvis)-0.51.*Karman-log(Karman);
end
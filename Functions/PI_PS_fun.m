function[PI_PS] = PI_PS_fun(p,n,freq,psi)

Karman = 0.41;

PI_PS=@(PI_PS_var) 2.*PI_PS_var-log(1+PI_PS_var)-Karman.*p.V./...
    p.u_tau-log(p.deltastar_PS.*p.V./p.kinvis)-0.51.*Karman-log(Karman);
end
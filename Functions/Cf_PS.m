function[Cf_SS] = Cf_PS(p,n,freq,psi)
    %p is the input parameters, and s is the section of the blade.

    Re_T = (p.density*p.U_e_PS(n)*p.c_R(n))/0.0000181206;
    Cf = 0.0986/(log10(Re_T)-1.22)^2;
end
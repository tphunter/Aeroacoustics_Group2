function[k_bar] = k_bar_fun(p,n,freq,psi)
    k_bar = freq.*p.chord(n)./(2.*p.c);
end


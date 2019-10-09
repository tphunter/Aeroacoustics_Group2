function[k_bar] = k_bar_fun(p,n)
    k_bar = @(freq) freq.*p.chord(n)./(2.*p.c);
end


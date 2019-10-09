function[mu_bar] = mu_bar_fun(p,n)
    beta_PG = sqrt(1-p.Mach^2); %Quizas cambiar
    mu_bar = @(freq) freq.*p.chord(n)./(2.*p.c.*beta_PG.^2);
end
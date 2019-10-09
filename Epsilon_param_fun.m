function[Epsilon_param] = Epsilon_param_fun(p,n)
kappa_bar = kappa_bar_fun(p,n);

Epsilon_param = @(freq) (exp(4.*1i.*kappa_bar).*(1-(1+1i).*fresnel(4.*kappa_bar)));
end
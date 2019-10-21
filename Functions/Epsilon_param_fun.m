function[Epsilon_param] = Epsilon_param_fun(p,n,freq,psi)
kappa_bar = kappa_bar_fun(p,n,freq,psi);

Epsilon_param =(exp(4.*1i.*kappa_bar).*(1-(1+1i).*fresnel(4.*kappa_bar)));
end
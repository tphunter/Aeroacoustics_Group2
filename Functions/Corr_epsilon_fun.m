function[Corr_epsilon] = Corr_epsilon_fun(p,n,freq,psi)
    epsilon = epsilon_fun(p,n,freq,psi);
    Epsilon_param = Epsilon_param_fun(p,n,freq,psi);
    
    Corr_epsilon = Epsilon_param+(epsilon-1)*imag(Epsilon_param);
end
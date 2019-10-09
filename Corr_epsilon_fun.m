function[Corr_epsilon] = Corr_epsilon_fun(p,n)
    epsilon = epsilon_fun(p,n);
    Epsilon_param = Epsilon_param_fun(p,n);
    
    Corr_epsilon = @(freq) Epsilon_param(n)+(epsilon(n)-1)*imag(Epsilon_param(n));
end
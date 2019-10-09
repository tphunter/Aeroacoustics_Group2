function[kappa_bar] = kappa_bar_fun(p,n)
    %p is the input parameters, and s is the section of the blade.
    mu_bar = mu_bar_fun(p,n);
    k_bar_y = k_bar_y_fun(p,n);
    beta_PG = sqrt(1-p.Mach.^2); %Quizas cambiar
    kappa_bar = @(freq) mu_bar.^2 - k_bar_y.^2./beta_PG.^2; %Check if needed @(freq)
end
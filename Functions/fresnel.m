function[E] = fresnel(x)
    fun = @(t) (exp(-1i.*t))./(sqrt(2.*pi.*t));
    E = integral(fun,0,x,'ArrayValued',1);
end
function[l_y] = l_y_fun(p,n)
    %p is the input parameters, and s is the section of the blade.
    
    l_y = @(freq) p.b*p.U_c/freq; %Maybe change M
end
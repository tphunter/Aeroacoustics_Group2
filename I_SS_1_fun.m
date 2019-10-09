function[I_SS_1] = I_SS_1_fun(p,n)
    %p is the input parameters, and s is the section of the blade.
    B_SS = B_SS_fun(p,n);
    C_SS = B_SS_fun(p,n);
    
    I_SS_1 = @(freq) 1i.*exp(2.*1i.*C_SS)./C_SS.*((1+1i).*...
        exp(-2.*1i.*C_SS).*sqrt(B_SS./(B_SS-C_SS)).*fresnel(2.*B_SS-2.*C_SS)...
        -(1+1i).*fresnel(2.*B_SS)+1);%k_x might be a vector; %Maybe change M
end
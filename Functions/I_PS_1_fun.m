function[I_PS_1] = I_PS_1_fun(p,n,freq,psi)
    %p is the input parameters, and s is the section of the blade.
    B_PS = B_PS_fun(p,n,freq,psi);
    C_PS = C_PS_fun(p,n,freq,psi);
    
    I_PS_1 = 1i.*exp(2.*1i.*C_PS)./C_PS.*((1+1i).*...
        exp(-2.*1i.*C_PS).*sqrt(B_PS./(B_PS-C_PS)).*fresnel(2.*B_PS-2.*C_PS)...
        -(1+1i).*fresnel(2.*B_PS)+1);%k_x might be a vector; %Maybe change M
end
function[phi_PP_PS] = phi_PP_PS_fun(p,n,freq,psi)
    
    PI_PS_function = PI_PS_fun(p,n,freq,psi);
    PI_PS = fzero(PI_PS_function,1);
    beta_c_PS = p.theta_PS(n)/p.tau_w_PS(n).*p.dP_dx_PS(n);
    Re_T_PS = p.u_tau.*p.delta_PS(n).*sqrt(p.C_f./2)./p.kinvis; %p.u_tau, p.C_f./2 might be array

    phi_PP_PS =(p.tau_w_PS(n).^2.*p.deltastar_PS(n)./(p.V).*0.78.*...
    (1.8.*PI_PS.*beta_c_PS+6).*(freq*p.deltastar_PS(n)./p.V).^2)./...
    (((freq.*p.deltastar_PS(n)./p.V).^0.75+0.105).^3.7+...
    (3.76.*Re_T_PS.^(-0.57).*(freq.*p.deltastar_PS(n)./p.V)).^7);
end
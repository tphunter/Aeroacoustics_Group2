function[phi_PP_SS] = phi_PP_SS_fun(p,n,freq,psi)
    
    PI_SS_function = PI_SS_fun(p,n,freq,psi);
    PI_SS = fzero(PI_SS_function,1);
    beta_c_SS = p.theta_SS(n)/p.tau_w_SS(n).*p.dP_dx_SS(n);
    Re_T_SS = p.u_tau.*p.delta_SS(n).*sqrt(p.C_f./2)./p.kinvis; %p.u_tau, p.C_f./2 might be array

    phi_PP_SS =(p.tau_w_SS(n).^2.*p.deltastar_SS(n)./(p.V).*0.78.*...
    (1.8.*PI_SS.*beta_c_SS+6).*(freq*p.deltastar_SS(n)./p.V).^2)./...
    (((freq.*p.deltastar_SS(n)./p.V).^0.75+0.105).^3.7+...
    (3.76.*Re_T_SS.^(-0.57).*(freq.*p.deltastar_SS(n)./p.V)).^7);
end
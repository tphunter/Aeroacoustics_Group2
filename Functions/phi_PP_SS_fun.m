function[phi_PP_SS] = phi_PP_SS_fun(p,n,freq,psi)
error=10;
PI_SS_var=-10;
PI_SS=10;
Karman=0.41;
Cf_SSvar = Cf_SS(p,n,freq,psi);
u_tau_SS = u_tau_SS(p,n,freq,psi);
while PI_SS>0.1
    PI_SS=2.*PI_SS_var-log(1+PI_SS_var)-Karman.*p.V./...
    p.u_tau+log(p.deltastar_SS.*p.V./p.kinvis)+0.51.*Karman+log(Karman);
    PI_SS_var=PI_SS_var+0.01;
   end
    beta_c_SS = p.theta_SS(n)/p.tau_w_SS(n).*p.dP_dx_SS(n);
    Re_T_SS = u_tau_SS.*p.delta_SS(n).*sqrt(Cf_SSvar./2)./p.kinvis; %p.u_tau, p.C_f./2 might be array
    phi_PP_SS =(p.tau_w_SS(n).^2.*p.deltastar_SS(n)./(p.V).*0.78.*...
    (1.8.*PI_SS.*beta_c_SS+6).*(freq*p.deltastar_SS(n)./p.V).^2)./...
    (((freq.*p.deltastar_SS(n)./p.V).^0.75+0.105).^3.7+...
    (3.76.*Re_T_SS.^(-0.57).*(freq.*p.deltastar_SS(n)./p.V)).^7);
end
function[phi_PP_PS] = phi_PP_PS_fun(p,n,freq,psi)
error=10;
PI_PS_var=100;
PI_PS=10;
Karman=0.41;
Cf_PSvar = Cf_PS(p,n,freq,psi);
u_tau_PSvar = u_tau_PS(p,n,freq,psi);
while abs(PI_PS)>10
    PI_PS=2.*PI_PS_var-log(1+PI_PS_var)-Karman.*p.V./...
    u_tau_PSvar+log(p.deltastar_PS(n).*p.V./p.kinvis)+0.51.*Karman+log(Karman);
    PI_PS_var=PI_PS_var+0.01;
   end
    beta_c_PS = p.theta_PS(n)/p.tau_w_PS(n).*p.dP_dx_PS(n);
    Re_T_PS = u_tau_PSvar.*p.delta_PS(n).*sqrt(Cf_PSvar./2)./p.kinvis; %u_tau_PSvar, p.C_f./2 might be array
    phi_PP_PS =(p.tau_w_PS(n).^2.*p.deltastar_PS(n)./(p.V).*0.78.*...
    (1.8.*PI_PS_var.*beta_c_PS+6).*(freq*p.deltastar_PS(n)./p.V).^2)./...
    (((freq.*p.deltastar_PS(n)./p.V).^0.75+0.105).^3.7+...
    (3.76.*Re_T_PS.^(-0.57).*(freq.*p.deltastar_PS(n)./p.V)).^7);
end
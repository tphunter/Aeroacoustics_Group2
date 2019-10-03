%% Parameters
%add all lengths in [m]
%all angles in [rad]
%all time in [s]
%all mass in [kg]
function p = Parameters()
%% POSITION
p.x = 10;
p.y = 10;
p.z = 10;
p.x_vector = [p.x,p.y,p.z];
p.i = [1,0,0];
p.j = [0,1,0];
p.k = [0,0,1];
p.R_0 = norm(x_vector);
p.theta =0 ;
p.gamma = asin(L/D);
p.phi = 0;



%% INPUT PARAMETERS
p.V=;
p.press=101325;
p.temp=;
p.density=p.press/(287*p.temp);
p.c=sqrt(p.temp*1.4*287);
p.Mach= p.V/p.c;
p.dynvis=;
p.kinvis=p.dynvis/p.density;
p.diam=;
p.R1=p.diam*0.5;
p.omega= (pi*6000)/30;
p.B=;
p.sections=;

p.r_R= [];
p.c_R= [];
p.chord=p.c_R*p.R1;
p.beta= [];
p.cl = [];
p.cd = [];
p.shapefactor= [];
p.U_e_SS= [];
p.dP_dx_SS= [];
p.delta_SS= [];
p.deltastar_SS= [];
p.theta_SS= [];
p.tau_w_SS= [];
p.U_e_PS= [];
p.dP_dx_PS= [];
p.delta_PS= [];
p.deltastar_PS= [];
p.theta_PS= [];
p.tau_w_PS= [];
p.k_x=p.omega/p.V;

%% TONAL



%% BROADBAND
p.U_c=;
p.u_tau=;
p.C_f=;
p.turb_kin_en=;
p.diss_rate=;
end
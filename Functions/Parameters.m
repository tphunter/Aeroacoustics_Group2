%% Parameters
%add all lengths in [m]
%all angles in [rad]
%all time in [s]
%all mass in [kg]
function p = Parameters()

fid = fopen('Data1.txt', 'rt');
data_cell = textscan(fid, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter',' ');
fclose(fid);

%% POSITION

p.R_0 = 1;
p.theta = pi/4;
p.cpcv=1.4;
p.phi = 0;
p.pref = 2E-5;
p.x = p.R_0*sin(p.phi)*cos(p.theta);
p.y = p.R_0*sin(p.phi)*sin(p.theta);
p.z = p.R_0*cos(p.phi);
p.x_vector = [p.x,p.y,p.z];
p.i = [1,0,0];
p.j = [0,1,0];
p.k = [0,0,1];




%% INPUT PARAMETERS
p.V= 18;
p.cpcv=1.4;
p.Rgas=287;
p.press=101325;
p.temp= 288.15;
p.density=p.press/(p.Rgas*p.temp);
p.c=sqrt(p.temp*p.cpcv*p.Rgas);
p.Mach= p.V/p.c;
p.dynvis=0.000018;
p.kinvis=p.dynvis/p.density;
p.diam=0.3556;
p.R1=p.diam*0.5;
p.omega= (pi*6000)/30;
p.B=2;
p.sections=30;


p.r_R = data_cell{1,1};
p.c_R = data_cell{1,2};
p.beta = data_cell{1,3};
p.cl = data_cell{1,4};
p.cd = data_cell{1,5};
p.shapefactor = data_cell{1,6};
p.U_e_SS = data_cell{1,7};
p.dP_dx_SS = data_cell{1,8};
p.delta_SS = data_cell{1,9};
p.deltastar_SS = data_cell{1,10};
p.theta_SS = data_cell{1,11};
p.tau_w_SS = data_cell{1,12};
p.U_e_PS = data_cell{1,13};
p.dP_dx_PS = data_cell{1,14};
p.delta_PS = data_cell{1,15};
p.deltastar_PS = data_cell{1,16};
p.theta_PS = data_cell{1,17};
p.tau_w_PS = data_cell{1,18};

p.gamma = asin(p.cl./p.cd);
p.chord=p.c_R*p.R1;
p.k_x=p.omega/p.V;

%% TONAL



%% BROADBAND
p.U_c=0.7*p.V;
p.u_tau=1;
p.C_f=1;
p.k_y = 1;
p.beta = 1;
p.b=1/2.1;
end
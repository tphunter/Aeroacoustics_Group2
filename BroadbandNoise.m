p = Parameters()



%% TRAILING EDGE

%1.- Radiation integral
k_bar=p.omega*p.chord/(2*a_0);
epsilon=1/sqrt(1+1/(4*k_bar));
sigma=sqrt(p.x^2+(1-p.Mach^2)*(p.y^2+p.z^2));
k_bar_y=k_y*p.chord/2;
k_bar_x=k_x*p.chord/2;
mu_bar=p.omega*p.chord/(2*p.c*beta^2); %This beta is not the angle I think 
B=k_bar_x-p.Mach*mu_bar+k_bar;
C=k_bar_x-mu_bar*(param.x/sigma-p.Mach);
theta=k_bar-mu_bar*param.x/sigma;
alpha=p.V/p.U_c;
H=((1+1i)exp(-4*1i*k_bar)(1-theta^2))/(2*sqrt(PI_SS)*(alpha-1)*k_x*sqrt(B));
G=(1+epsilon)exp(1i(2*k_bar+theta))sin(theta-2*k_bar)/(theta-2*k_bar)+(1-epsilon)*exp(1i(-2*k_bar+theta))sin(theta+2*k_bar)/(theta+2*k_bar)+(1+epsilon)(1-1i)/(2*(theta-2*k_bar))exp(4*1i*k_bar)*fresnelc(4*k_bar)-(1-epsilon)(1+1i)/(2*(theta+2*k_bar))exp(4*1i*k_bar)*fresnelc(4*k_bar)+exp(2*1i*theta)/2*sqrt(2*k_bar/theta)*fresnelc(2*theta)((1-epsilon)(1+1i)/(theta+2*k_bar)-(1+epsilon)(1+1i)/(theta-2*k_bar)); %I use Fresnel cosine func
I_TE_1=1i*e^(2*1i*C)/C*((1+1i)*exp(-2*1i*C)*sqrt(B/(B-C))*fresnel(2*B-2*C)-(1+1i)*fresnel(2*B)+1); %I use the Fresnel integral of cosine, ask if we should use the sine one.
I_TE_2=H*(exp(4*1i*k_bar)(1-(1+1i)*fresnel(4*k_bar)))^epsilon+H(exp(-2*1i*theta)+1i*(theta+k_bar_x+p.Mach*mu_bar-k_bar)*G); %Take a look at definition of epsilon (slide 8 lecture 4), it should not be an exponent.

I_TE=I_TE_1+I_TE_2;

%2.- Wall-pressure wave-number frequency spectrum
b = 1/2.1; %Calibration constant for l_y
kappa=0.41;

%SS
Re_T_SS=p.u_tau*p.delta_SS*sqrt(p.C_f/2)/p.kinvis;
beta_c_SS=theta/p.tau_w_SS*p.dP_dx_SS;
fun_PI_SS=@(PI_SS_var) 2*PI_SS_var-log(1+PI_SS_var)-kappa*p.V/p.u_tau-log(p.deltastar_SS*p.V/p.kinvis)-0.51*kappa-log(kappa);
PI_SS=fzero(fun_PI_SS,1);
phi_pp_omega_SS=(p.tau_w_SS^2*p.delta_star_SS/(p.V)0.78(1.8*PI_SS*beta_c_SS+6)(p.omega*p.delta_star_SS/p.V)^2)/(((p.omega*p.delta_star_SS/p.V)^0.75+0.105)^3.7+(3.76*Re_T_SS^(-0.57)(p.omega*p.delta_star_SS/p.V))^7);
l_y_TE_SS=b*p.U_c/p.omega;
phi_pp_SS=phi_pp_omega_SS*l_y_TE/PI_SS;
S_pp_TE_SS=(k_bar*p.z/(2*PI_SS*sigma^2))^2*2*p.chord*integral(phi_pp_SS*sin(p.R1/p.chord*(k_bar_y-k_bar*p.z/sigma))^2/((k_bar_y-k_bar*p.z/sigma)^2)*abs(I_TE)^2,-inf,inf);

%PS
Re_T_PS=p.u_tau*p.delta_PS*sqrt(p.C_f/2)/p.kinvis;
beta_c_PS=theta/p.tau_w_PS*p.dP_dx_PS;
fun_PI_PS=@(PI_PS_var) 2*PI_PS_var-log(1+PI_PS_var)-kappa*p.V/p.u_tau-log(p.deltastar_PS*p.V/p.kinvis)-0.51*kappa-log(kappa);
PI_PS=fzero(fun_PI_PS,1);
phi_pp_omega_PS=(p.tau_w_PS^2*p.delta_star_PS/(p.V)0.78(1.8*PI_PS*beta_c_PS+6)(p.omega*p.delta_star_PS/p.V)^2)/(((p.omega*p.delta_star_PS/p.V)^0.75+0.105)^3.7+(3.76*Re_T_PS^(-0.57)(p.omega*p.delta_star_PS/p.V))^7);
l_y_TE_PS=b*p.U_c/p.omega;
phi_pp_PS=phi_pp_omega_PS*l_y_TE/PI_PS;
S_pp_TE_PS=(k_bar*p.z/(2*PI_PS*sigma^2))^2*2*p.chord*integral(phi_pp_PS*sin(p.R1/p.chord*(k_bar_y-k_bar*p.z/sigma))^2/((k_bar_y-k_bar*p.z/sigma)^2)*abs(I_TE)^2,-inf,inf);
%% LEADING EDGE

%1.- Radiation integral
theta_plus=k_bar+mu_bar/sigma;
theta_minus=k_bar-mu_bar/sigma;
theta_bar=k_bar_x*(1-p.Mach*p.x/sigma)/(1-p.Mach^2)-pi/4;
I_LE_1=-1/pi*sqrt(2/(theta_minus*(k_bar_x+(1-p.Mach^2)*k_bar*theta_minus)))*exp(-1i*theta)*fresnel(2*theta);
I_LE_2=exp(-1i*theta_bar)/(pi*theta*sqrt(2*pi*(k_bar_x+(1-p.Mach^2)k_bar)))(1i*(1-exp(2*1i*theta_minus))-(1+1i)*(fresnel(4*k_bar)-exp(2*1i*theta_minus)*sqrt(2*k_bar/theta_plus)*fresnel(2*theta_plus)));

I_LE=I_LE_1+I_LE_2;

%2.- Phi omega omega
L_t=0.4*p.turb_kin_en^1.5/p.diss_rate;
u_bar_prime=sqrt(2*p.turb_kin_en/3);
k_e=sqrt(pi)/L_t*gamma(5/6)/gamma(1/3);
phi_omegaomega_omega = u_bar_prime^2*L_t/(2*pi*p.V)(1+8/3(p.k_x+k_e)^2)/((1+(p.k_x/k_e)^2)^(211/6));
l_y_LE=(gamma(1/3)/gamma(5/6))^2*(p.k_x/k_e)^2/((3+8*(p.k_x+k_e)^2)*sqrt(1+(p.k_x+k_e)^2));

phi_omegaomega=p.V*phi_omegaomega_omega*l_y_LE/pi;

%3.- Power spectral density
S_pp_LE=(p.density*k_bar*p.z/(sigma^2))^2*pi*p.V*p.R1/2*integral(phi_omegaomega*sin(p.R1/2*(k_y/sigma-k_y))^2/(PI_SS/2*(k_y/sigma-k_y)^2*abs(I_LE)^2),-inf,inf);


%% SUBFUNCTIONS
function[E] = fresnel(x)

fun = @(t) (exp(-1i*t))/(sqrt(2*PI_SS*t));
E = integral(fun,0,x);

end
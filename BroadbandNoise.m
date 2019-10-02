p = Parameters()



%% TRAILING EDGE

%1.- Radiation integral
k_bar=p.omega*p.chord/(2*a_0);
epsilon=1/sqrt(1+1/(4*k_bar));
sigma=sqrt(p.x^2+(1-p.Mach^2)*(p.y^2+p.z^2));
k_y_bar=k_y*p.chord/2;
k_x_bar=k_x*p.chord/2;
Avg_Mach=p.V/p.a_0;
mu_bar=p.omega*p.chord/(2*p.a_0*beta^2); %This beta is not the angle I think 
B=k_bar_x-Avg_Mach*mu_bar+k_bar;
C=k_bar_x-mu_bar*(p.x/sigma-Avg_Mach);
theta=k_bar-mu_bar*p.x/sigma;
H=((1+i)*exp(-4*i*k_bar)*(1-theta^2))/(2*sqrt(pi)*(alpha-1)*k_x*sqrt(B));
G=(1+epsilon)*exp(i*(2*k_bar+theta))*sin(theta-2*k_bar)/(theta-2*k_bar)+(1-epsilon)*exp(i*(-2*k_bar+theta))*sin(theta+2*k_bar)/(theta+2*k_bar)+(1+epsilon)*(1-i)/(2*(theta-2*k_bar))*exp(4*i*k_bar)*fresnelc(4*k_bar)-(1-epsilon)*(1+i)/(2*(theta+2*k_bar))*exp(4*i*k_bar)*fresnelc(4*k_bar)+exp(2*i*theta)/2*sqrt(2*k_bar/theta)*fresnelc(2*theta)*((1-epsilon)*(1+i)/(theta+2*k_bar)-(1+epsilon)*(1+i)/(theta-2*k_bar)); %I use Fresnel cosine func
I_TE_1=i*e^(2*i*C)/C*((1+i)*exp(-2*i*C)*sqrt(B/(B-C))*fresnelc(2*B-2*C)-(1+i)*fresnelc(2*B)+1); %I use the Fresnel integral of cosine, ask if we should use the sine one.
I_TE_2=H*(exp(4*i*k_bar)*(1-(1+i)*fresnelc(4*k_bar)))^epsilon+H*(exp(-2*i*theta)+i*(theta+k_bar_x+Avg_Mach*mu_bar-k_bar)*G); %Take a look at definition of epsilon (slide 8 lecture 4), it should not be an exponent.

I_TE=I_TE_1+I_TE_2;

%2.- Wall-pressure wave-number frequency spectrum
b=; %Calibration constant for l_y
kappa=0.41;
%Need to compute nu and C_f
Re_T=u_tau*p.delta_star*sqrt(C_f/2)/nu;
beta_c=theta/p.tau_w*p.dP_dx;
fun_PI=@(PI_var) 2*PI_var-log(1+PI_var)-kappa*p.V/u_tau-log(delta_star*p.V/nu)-0.51*kappa-log(kappa);
PI=fzero(fun_PI,1);
phi_pp_omega=(p.tau_w^2*p.delta_star/(p.V)*0.78*(1.8*PI*beta_c+6)*(p.omega*p.delta_star/p.V)^2)/(((p.omega*p.delta_star/p.V)^0.75+0.105)^3.7+(3.76*Re_T^(-0.57)*(p.omega*p.delta_star/p.V))^7);
l_y_TE=b*U_c/p.omega;
phi_pp=phi_pp_omega*l_y_TE/pi;

%3.- Power spectral density
S_pp_TE=(k_bar*p.z/(2*pi*sigma^2))^2*2*p.chord*integral(phi_pp*sin(param_L/param_chord*(k_bar_y-k_bar*p.z/sigma))^2/((k_bar_y-k_bar*p.z/sigma)^2)*abs(I_TE)^2,-inf,inf);

%% LEADING EDGE

%1.- Radiation integral
theta_plus=k_bar+mu_bar/sigma;
theta_minus=k_bar-mu_bar/sigma;
theta_bar=k_bar_x*(1-Avg_Mach*p.x/sigma)/(beta^2)-pi/4;
I_LE_1=-1/pi*sqrt(2/(theta_minus*(k_bar_x+beta^2*k_bar)))*exp(-i*theta)*fresnelc(2*theta);
I_LE_2=exp(-i*theta_bar)/(pi*theta*sqrt(2*pi*(k_bar_x+beta^2*k_bar)))*(i*(1-exp(2*i*theta_minus))-(1+i)*(fresnelc(4*k_bar)-exp(2*i*theta_minus)*sqrt(2*k_bar/theta_plus)*fresnelc(2*theta_plus)));

I_LE=I_LE_1+I_LE_2;

%2.- Phi omega omega
L_t=0.4*turb_kin_en^1.5/diss_rate;
u_bar_prime=sqrt(2*turb_kin_en/3);
k_e=sqrt(pi)/L_t*tau(5/6)/tau(1/3); %I DONT KNOW WHAT IS TAU, MIGHT BE SOME FUNCTION
phi_omegaomega_omega=u_bar_prime^2*L_t/(2*pi*p.V)*(1+8/3*(k_x+k_e)^2)/((1+(k_x/k_e)^2)^(211/6))
l_y_LE=(tau(1/3)/tau(5/6))^2*(k_x/k_e)^2/((3+8*(k_x+k_e)^2)*sqrt(1+(k_x+k_e)^2));

phi_omegaomega=p.V*phi_omegaomega_omega*l_y_LE/pi;

%3.- Power spectral density
S_pp_LE=(p.rho_0*k_bar*p.z/(sigma^2))^2*pi*p.V*p.L/2*integral(phi_omegaomega*sin(p.L/2*(k*p.y/sigma-k*p.y))^2/(pi/2*(k*p.y/sigma-k*p.y)^2*abs(I_LE)^2),-inf,inf);
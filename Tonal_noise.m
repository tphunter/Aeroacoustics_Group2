%% Tonal Noise Functions
p = Parameters()
sres = p.sres;
phi = p.phi ;
Omega = p.Omega;

s=50;
P_mB= i*m*B^2*Omega/(4*pi*c*R_0)*sum_-s^s*F_s*exp(-i*(m*B-s)*pi/2)*exp(i*(m*B-s)(phi-Omega_s*R_0/c))*JmB-s(m*B*M*sin(theta))(-(m*B-s)/(m*B)*sin(gamma)/M+cos(theta)*cos(gamma));
Omegas = m*B*Omega/(m*B-s);

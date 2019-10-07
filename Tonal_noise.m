%% Tonal Noise Functions
p = Parameters();
a = Get_P_mB(p);

function P_mB = Get_P_mB(p)
Omega = p.omega;
theta = p.theta;
phi = p.phi;
R_0 = p.R_0;
c = p.c;
m=1;
B=p.B;
M = p.R1*Omega/p.c;
mB=m*B;
cl = p.cl;
cd = p.cd;
force = sqrt(cd.^2+cl.^2);
gamma = atan(cd./cl);


F_s = fft(force);

P_mB = [];
s=0.5*length(F_s);

for si=-s:1:s
    Omega_s = mB*Omega/(mB-si);
    F_si = F_s(0.5*length(F_s)+si+1);
    P_mBi = F_si...
        *exp(-1i*(mB-si)*pi/2)*exp(1i*(mB-si)*(phi-Omega_s*R_0/c))...
        *besselj(mB-si,mB*M*sin(theta))*...
        (-(mB-si)/(mB)*sin(gamma)/M+cos(theta)*cos(gamma));
    P_mBi
    p_mB = P_mB + P_mBi;
end
P_mB = 1i*mB*B*Omega/(4*pi*c*R_0)*p_mB;
end


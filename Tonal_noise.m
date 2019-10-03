%% Tonal Noise Functions
p = Parameters();
a = Get_P_mB(p);

function P_mB = Get_P_mB(p)
Omega = p.Omega;
theta = p.theta;
gamma = p.gamma;
phi = p.phi;
R_0 = p.R_0;
c = p.c;
m=p.m;
B=p.B;
M = p.M;
mB=m*B;
F_s = fft(p.force);

P_mB = [];
s=50;

for s=-s:1:s
    Omega_s = mB*Omega/(mB-s);
    F_si = F_s(s);
    P_mBi = F_si...
        *exp(-1i*(mB-s)*pi/2)*exp(1i*(mB-s)*(phi-Omega_s*R_0/c))...
        *jbessel(mB-s,mB*M*sin(theta))*...
        (-(mB-s)/(mB)*sin(gamma)/M+cos(theta)*cos(gamma));
    p_mB = P_mB + P_mBi;
end
P_mB = 1i*mB*B*Omega/(4*pi*c*R_0)*p_mB;
end


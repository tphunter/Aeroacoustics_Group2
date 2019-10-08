%% Tonal Noise Functions
p = Parameters();
a = Get_P_mB(p);
freq = real(a);
ampl = imag(a);
plot( freq, ampl, 'bo')


function P_mBfinal = Get_P_mB(p)
Omega = p.omega;
theta = p.theta;
phi = p.phi;
R_0 = p.R_0;
c = p.c;
m=1;
B=p.B;

mB=m*B;
cl = p.cl;
cd = p.cd;
force = sqrt(cd.^2+cl.^2);
gamma = atan(cd./cl);


P_mBfinal = [];
for i = 1:length(force)
    P_mB = 0;
    force_i = force(i) + 0.1 *(0.5-rand(10,1));
    F_s = fft(force_i,501);

    M = p.r_R(i)*p.R1*Omega/p.c;
    s=floor(0.5*length(F_s));
    p_mB = 0;
    for si=-s:1:s
        if si ~= mB
        Omega_s = mB*Omega/(mB-si);
        F_si = F_s(s+si+1);
        gamma_i = gamma(i);
        P_mBi = F_si...
            *exp(-1i*(mB-si)*pi/2)*exp(1i*(mB-si)*(phi-Omega_s*R_0/c))...
            *besselj(mB-si,mB*M*sin(theta))*...
            (-(mB-si)/(mB)*sin(gamma_i)/M+cos(theta)*cos(gamma_i));
      
        p_mB = p_mB+ P_mBi;
        
        else
            si = mB + 10^-9;
            Omega_s = mB*Omega/10^-9;
        F_si = F_s(s+mB+1);
        gamma_i = gamma(i);
        P_mBi = F_si...
            *exp(-1i*(mB-si)*pi/2)*exp(1i*(mB-si)*(phi-Omega_s*R_0/c))...
            *besselj(mB-si,mB*M*sin(theta))*...
            (-(mB-si)/(mB)*sin(gamma_i)/M+cos(theta)*cos(gamma_i));
        
        p_mB = p_mB+ P_mBi;
        end
    end
    P_mB = 1i*mB*B*Omega/(4*pi*c*R_0)*p_mB;
    P_mBfinal = [P_mBfinal,P_mB];
end
end

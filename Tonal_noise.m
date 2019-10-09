%% Tonal Noise Functions
p = Parameters();
% m = 2000;
% a = Get_P_mB(p,m);
% freq = abs(real(a));
% ampl = abs(imag(a));
% %freq = real(a);
% %ampl = imag(a);
% plot( freq, ampl , 'x')
P=[];
mobj = 500;
for m = 1:mobj
    P = [P ,Get_P_mB(p,m)];
end
      

        
plot(linspace(0,mobj*p.B,length(P)), abs(P),'x')

function P_mBfinal = Get_P_mB(p,m)
Omega = p.omega;
theta = p.theta;
phi = p.phi;
R_0 = p.R_0;
c = p.c;
%m=1;
B=p.B;

mB=m*B;
cl = p.cl;
cd = p.cd;
force = sqrt(cd.^2+cl.^2);
gamma = atan(cd./cl);
resolution = 101;
n= resolution;
t = linspace(0,2*pi/Omega,n);

P_mBfinal = [];
for i = 1:length(force)
    P_mB = 0;
    F_s = 0;
    force_i = 1 + 1*p.r_R(i)*(normrnd(0,1,[n,1]));%*sin(t*4*Omega)
%     plot(t,force_i)
    
    F_s = fft(force_i,n);
    f = (0:length(F_s)-1)*300/length(F_s);
%     plot(f, abs(F_s))
    F_s = abs(fftshift(F_s));
    fshift = (-n/2:n/2-1)*(300/n);    
%     plot(fshift, F_s)
    
    gamma_i = gamma(i);
    M = p.r_R(i)*p.R1*Omega/p.c;
    s=floor(0.5*length(fshift));
    p_mB = 0;
    for si=-s:1:s
        if si ~= mB
        Omega_s = mB*Omega/(mB-si);
        F_si = F_s(s+si+1);
        
        P_mBi = F_si...
            *exp(-1i*(mB-si)*pi/2)*exp(1i*(mB-si)*(phi-Omega_s*R_0/c))...
            *besselj(mB-si,mB*M*sin(theta))*...
            (-(mB-si)/(mB)*sin(gamma_i)/M+cos(theta)*cos(gamma_i));
      
        p_mB = p_mB+ P_mBi;
        end
        if si == mB
            si = mB + 10^-9;
            Omega_s = mB*Omega/10^-9;
        F_si = F_s(s+mB+1);
        
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

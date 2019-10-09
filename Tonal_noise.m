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
for m = 1:200
    P = [P ,Get_P_mB(p,m)];
end
      

        
plot(abs(real(P)), abs(imag(P)),'x')

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
resolution = 201;
n= resolution;
t = 0:4*pi/Omega/n:4*pi/Omega-4*pi/Omega/n;

P_mBfinal = [];
for i = 1:length(force)
    P_mB = 0;
    F_s = 0;
    force_i = force(i) + 0.1*sin(2*pi*t*30);%*p.r_R(i); %*(0.5-rand(n,1));%
    plot(linspace(0,2*pi/Omega-2*pi/Omega/n,length(force_i)),force_i)
    fshift = (-n/2:n/2-1)*(400/n);
    F_s = fft(force_i,n);
%     plot(fshift, F_s)
    F_s = abs(fftshift(F_s));    
    plot(fshift, F_s)
    
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

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
mobj = 1;
for m = 1:mobj
    P = [P ,Get_P_mB(p,m)];
end
      

        
plot(linspace(0,mobj*p.B*p.omega/2/pi,length(P)),abs(P))% 20*log10(abs(P)/20/10^-6))

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

n= 1001;
Fs = 2.1*5000;
T = 1/Fs;
L = n;
t = (0:L-1)*T;%linspace(0,10*2*pi/Omega,n);

P_mBfinal = [];
for i = 1:length(force)
    P_mB = 0;
    F_s = 0;
    force_i = 0.5*p.density*p.V^2*force(i)*p.c_R(i)^2/p.shapefactor(i);
    force_i = force_i - 0.5*sin(2*pi*75*t);%1.5*(1-p.r_R(i))*sin(t*2*pi*100).^2;%(normrnd(0,6,[n,1]));%*sin(t*4*Omega)1*p.r_R(i)*
%     plot(t(1:500),force_i(1:500))
    
    F_s = fft(force_i,n);
    f = (0:n-1)*(Fs/n);
%     plot(f, abs(F_s/(2*pi/Omega)))
    F_s = abs(fftshift(F_s));
    fshift = (-n/2:n/2-1)*(Fs/n);
%     plot(fshift, abs(F_s))
    
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
            *besselj(mB-si,mB*M*sin(theta),1)*...
            (-(mB-si)/(mB)*sin(gamma_i)/M+cos(theta)*cos(gamma_i));
      
        p_mB = p_mB+ P_mBi;
        end
        if si == mB
            
            Omega_s = Omega;
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

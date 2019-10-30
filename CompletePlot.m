clc;
clear;
warning('off','all');
%% Tonal Noise Functions
p = Parameters();
% m = 2000;
% a = Get_P_mB(p,m);
% freq = abs(real(a));
% ampl = abs(imag(a));
% %freq = real(a);
% %ampl = imag(a);
% plot( freq, ampl , 'x')

Pdir=[];
P=[];
P2 = [];
mobj = 27;

for m = 1:mobj
    P = [P ,sum(Get_P_mB(p,m, p.theta, p.phi,0.1))];
    P2 = [P2 ,sum(Get_P_mB(p,m, p.theta, p.phi,0))];
end  
% mobj = 27;
% P2 = [];
% for m = 1:mobj
%     P2 = [P2 ,sum(Get_P_mB(p,m, p.theta, p.phi,0))];
% end

lengthpsi=4;

i_psi = 2; 
i_freq = 1;
i_sect=2;
for freqnondop=100:1000:5100
    tic
    for n=1:p.sections-1
        for psi=0:2*pi/lengthpsi:2*pi
            
            Mt = p.omega*p.R1/p.c ;
            Mz = p.Mach;
            freqratio = 1 + Mt*sin(p.theta)*sin(psi)/(sqrt(1-Mz^2*sin(p.theta)^2));
            freq=freqratio*freqnondop;
            Spp_SS(i_psi,i_freq,i_sect-1) = Spp_SS_fun(p,n,freq,psi);
            Spp_PS(i_psi,i_freq,i_sect-1) = Spp_PS_fun(p,n,freq,psi);
            
            i_psi = i_psi + 1;
        end
        for i_psi=1:lengthpsi-1
            fun_psi_SS(i_psi) = 2*pi/lengthpsi*(freqratio*Spp_SS(i_psi,i_freq,i_sect-1)+freqratio*Spp_SS(i_psi+1,i_freq,i_sect-1));
            fun_psi_PS(i_psi) = 2*pi/lengthpsi*(freqratio*Spp_PS(i_psi,i_freq,i_sect-1)+freqratio*Spp_PS(i_psi+1,i_freq,i_sect-1));
        end    
        Spp_SS_sect_freq(i_freq,i_sect-1) = p.B/(2*pi)*sum(fun_psi_SS);
        Spp_PS_sect_freq(i_freq,i_sect-1) = p.B/(2*pi)*sum(fun_psi_PS);
            
        Spp_sect_freq(i_freq,i_sect-1)=Spp_SS_sect_freq(i_freq,i_sect-1)+Spp_PS_sect_freq(i_freq,i_sect-1);
        i_sect=i_sect+1;
    end
    Spp_freq_rad(i_freq)=0.259/2*(Spp_sect_freq(i_freq,i_sect-2)+Spp_sect_freq(i_freq,i_sect-3));
    i_freq=i_freq+1;
    toc
end

freqSpp=100:1000:5100;
Spp_freq=1/(2*pi)*Spp_freq_rad;
pol_Spp_freq=polyfit(freqSpp,Spp_freq,length(freqSpp)-1);
pol_Spp_freq_fun=@(m) pol_Spp_freq(1).*m.^5+pol_Spp_freq(2).*m.^4+pol_Spp_freq(3).*m.^3+pol_Spp_freq(4).*m.^2+pol_Spp_freq(5).*m.^2+pol_Spp_freq(6);

index=1;
for freqSPL=100:10:4990
    SPL(index)=20*log10(sqrt(abs(integral(pol_Spp_freq_fun,freqSPL,freqSPL+10)))/(2E-5));
    index=index+1;
end
freqSPL=100:10:4990;

figure(1)
stem(linspace(200,mobj*200,length(P)),20*log10(abs(P2)/2/10^-5))
hold on
plot(freqSPL,SPL)
title('Noise Power vs Frequency')
xlabel('Frequency')
ylabel('Noise Power [dB]')


figure(2)
stem(linspace(200,mobj*200,length(P)),20*log10(abs(P)/2/10^-5))
hold on
plot(freqSPL,SPL)
title('Noise Power vs Frequency')
xlabel('Frequency')
ylabel('Noise Power [dB]')



function P_mBfinal = Get_P_mB(p,m, theta, phi,k)
Omega = p.omega;
%theta = p.theta;
%phi = p.phi;
R_0 = p.R_0;
c = p.c;
%m=1;
B=p.B;

mB=m*B;
cl = p.cl;
cd = p.cd;
force = sqrt(cd.^2+cl.^2);
gamma = atan(cd./cl);

n= 10001;
Fs = 2.1*5000;
T = 1/Fs;
L = n;
t = (0:L)*T;%linspace(0,10*2*pi/Omega,n);

P_mBfinal = [];
for i = 1:length(force) %for i-th panel
    P_mB = 0;
    F_s = 0;
    force_i = force(i);%0.5*p.density*p.V^2*force(i)*p.c_R(i)^2/p.shapefactor(i); %change V, 
    force_i = force_i - k*sin(2*pi*170*t);%1.5*(1-p.r_R(i))*sin(t*2*pi*100).^2;%(normrnd(0,6,[n,1]));%*sin(t*4*Omega)1*p.r_R(i)*
%     plot(t(1:500),force_i(1:500))
    
    F_s = fft(force_i,n);
    f = (0:n-1)*(Fs/n);
%     plot(f, abs(F_s/(2*pi/Omega)))
    F_s = abs(fftshift(F_s))/n;
    fshift = (-n/2:n/2-1)*(Fs/n);
%     plot(linspace(-5250,5250,n), abs(F_s))
    
    gamma_i = gamma(i);
    M = p.r_R(i)*p.R1*Omega/p.c;
    s=floor(0.5*n);
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
        F_si = F_s(s+si+1);
        
        P_mBi = F_si...
            *exp(-1i*(mB-si)*pi/2)*exp(1i*(mB-si)*(phi-Omega_s*R_0/c))...
            *besselj(mB-si,mB*M*sin(theta),1)*...
            (-(mB-si)/(mB)*sin(gamma_i)/M+cos(theta)*cos(gamma_i));
        
        p_mB = p_mB+ P_mBi;
        end
    end
    P_mB = 1i*mB*B*Omega/(4*pi*c*R_0)*p_mB;
    P_mBfinal = [P_mBfinal,P_mB];
end
end

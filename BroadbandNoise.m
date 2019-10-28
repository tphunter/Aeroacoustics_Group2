clc;
clear;
warning('off','all');
p = Parameters();
i_psi = 2;
i_freq = 1;
i_sect=1;
lengthpsi=4;

for freqnondop=100:1000:5100
    tic
    for psi=0:2*pi/lengthpsi:2*pi
        for n=1:p.sections-1
            
            Mt = p.omega*p.r_R(n)*p.R1/p.c ;
            Mz = p.Mach;
            freqratio = 1 + Mt.*sin(p.theta).*sin(psi)/(sqrt(1-Mz.^2*sin(p.theta).^2));
            freq=freqratio.*freqnondop;
            Spp_SS_sect(i_psi,i_freq,i_sect) = Spp_SS_fun(p,n,freq,psi);
            Spp_PS_sect(i_psi,i_freq,i_sect) = Spp_PS_fun(p,n,freq,psi);
            i_sect=i_sect+1;
            
        end
        Spp_SS(i_psi,i_freq)=0;
        Spp_PS(i_psi,i_freq)=0;
        for i_sect=1:length(n)-1 
            Spp_SS(i_psi,i_freq)=Spp_SS(i_psi,i_freq)+Spp_SS_sect(i_psi,i_freq,i_sect);
            Spp_PS(i_psi,i_freq)=Spp_PS(i_psi,i_freq)+Spp_PS_sect(i_psi,i_freq,i_sect);
        end
        i_psi = i_psi + 1;
    end
    for i_psi=1:length(psi)-1
        fun_psi_SS(i_psi) = 2*pi/lengthpsi*freqratio*(Spp_SS(i_psi,i_freq)+Spp_SS(i_psi+1,i_freq));
        fun_psi_PS(i_psi) = 2*pi/lengthpsi*freqratio*(Spp_PS(i_psi,i_freq)+Spp_PS(i_psi+1,i_freq));
    end
    Spp_SS_sect_freq(i_freq,i_sect) = p.B/(2*pi)*sum(fun_psi_SS);
    Spp_PS_sect_freq(i_freq,i_sect) = p.B/(2*pi)*sum(fun_psi_PS);
    
    Spp_sect_freq(i_freq,i_sect)=Spp_SS_sect_freq(i_freq,i_sect)+Spp_PS_sect_freq(i_freq,i_sect-1);
    Spp_freq_rad(i_freq)=0.259/2*(Spp_sect_freq(i_freq,i_sect)+Spp_sect_freq(i_freq,i_sect+1));
    i_freq=i_freq+1;
    toc
end

freqSpp=100:1000:5100;
Spp_freq=1/(2*pi)*Spp_freq_rad;
pol_Spp_freq=polyfit(freqSpp,Spp_freq,length(freqSpp)-1);
pol_Spp_freq_fun=@(m) pol_Spp_freq(1).*m.^5+pol_Spp_freq(2).*m.^4+pol_Spp_freq(3).*m.^3+pol_Spp_freq(4).*m.^2+pol_Spp_freq(5).*m.^2+pol_Spp_freq(6);

index=1;
for freqSPL=100:10:4990
    SPL(index)=20*log10(sqrt(integral(pol_Spp_freq_fun,freqSPL,freqSPL+10))/(2E-5));
    index=index+1;
end
freqSPL=100:10:4990;

figure(1)
plot(freqSPL,SPL)
xlabel('Frequency (Hz)')
ylabel('SPL (dB)')

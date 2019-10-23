clc;
clear;
warning('off','all');
p = Parameters();
i_psi = 2;
i_freq = 1;
i_sect=2;
lengthpsi=4;


for freq=0.01:5000:15000.01
    tic
    for n=1:p.sections-1
        tic
        for psi=0:2*pi/lengthpsi:2*pi
            
            Mt = p.omega*p.R1/p.c ;
            Mz = p.Mach;
            Thetha(i_psi,i_freq,i_sect-1) = Theta_minus_fun(p,n,freq,psi);
            freqratio(i_psi,i_freq,i_sect-1) = 1 + Mt*sin(Thetha(i_psi,i_freq,i_sect-1))*sin(psi)/(sqrt(1-Mz^2*sin(Thetha(i_psi,i_freq,i_sect-1))^2));
            freqratio(i_psi-1,i_freq,i_sect-1) = 1 + Mt*sin(Thetha(i_psi-1,i_freq,i_sect-1))*sin(psi)/(sqrt(1-Mz^2*sin(Thetha(i_psi-1,i_freq,i_sect-1))^2));
            Spp_SS(i_psi,i_freq,i_sect-1) = Spp_SS_fun(p,n,freq,psi);
            Spp_PS(i_psi,i_freq,i_sect-1) = Spp_PS_fun(p,n,freq,psi);
            
            i_psi = i_psi + 1;
        end
        toc
        for i_psi=1:lengthpsi-1
            fun_psi_SS(i_psi) = 2*pi/lengthpsi*(freqratio(i_psi,i_freq,i_sect-1)*Spp_SS(i_psi,i_freq,i_sect-1)+freqratio(i_psi+1,i_freq,i_sect-1)*Spp_SS(i_psi+1,i_freq,i_sect-1));
            fun_psi_PS(i_psi) = 2*pi/lengthpsi*(freqratio(i_psi,i_freq,i_sect-1)*Spp_PS(i_psi,i_freq,i_sect-1)+freqratio(i_psi+1,i_freq,i_sect-1)*Spp_PS(i_psi+1,i_freq,i_sect-1));
        end    
        Spp_SS_sect_freq(i_freq,i_sect-1) = p.B/(2*pi)*sum(fun_psi_SS);
        Spp_PS_sect_freq(i_freq,i_sect-1) = p.B/(2*pi)*sum(fun_psi_PS);
            
        Spp_sect_freq(i_freq,i_sect-1)=Spp_SS_sect_freq(i_freq,i_sect-1)+Spp_PS_sect_freq(i_freq,i_sect-1);
        i_sect=i_sect+1;
    end
    Spp_freq(i_freq)=0.259/2*(Spp_sect_freq(i_freq,i_sect-2)+Spp_sect_freq(i_freq,i_sect-3));
    i_freq=i_freq+1;
    toc
end

freq=100:1000:5000;
plot(freq,Spp_freq)


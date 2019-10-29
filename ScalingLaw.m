clc;
clear;
warning('off','all');
addpath('Functions')
p = Parameters();
lengthpsi=4;
theta = linspace(0, 2*pi,15);

R_0 = p.R_0;
n_iter = 1;

for R = R_0:25:100
    p.R_0 = R;
    for i=1:length(theta)
    i_psi = 2; 
    i_freq = 1;
    i_sect=2;
    p.theta=theta(i);
    for freq=100:1000:5100
        tic
        for n=1:p.sections-1
            for psi=0:2*pi/lengthpsi:2*pi

                Mt = p.omega*p.R1/p.c ;
                Mz = p.Mach;
                freqratio = 1 + Mt*sin(p.theta)*sin(psi)/(sqrt(1-Mz^2*sin(p.theta)^2));
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

    freqSPL=100:10:4990;

    OASPL=20*log10(sqrt(abs(integral(pol_Spp_freq_fun,100,5000)))/(2E-5));
    end
    RList(n_iter) = R;
    OASPLList(n_iter) = OASPL;
    n_iter = n_iter + 1;
    
end

Diam_0 = p.diam;
n_iter = 1;
for Diam = Diam_0:5:20
    p.diam = Diam;
    for i=1:length(theta)
    i_psi = 2; 
    i_freq = 1;
    i_sect=2;
    p.theta=theta(i);
    for freq=100:1000:5100
        tic
        for n=1:p.sections-1
            for psi=0:2*pi/lengthpsi:2*pi

                Mt = p.omega*p.R1/p.c ;
                Mz = p.Mach;
                freqratio = 1 + Mt*sin(p.theta)*sin(psi)/(sqrt(1-Mz^2*sin(p.theta)^2));
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


    freqSPL=100:10:4990;

    OASPL=20*log10(sqrt(abs(integral(pol_Spp_freq_fun,100,5000)))/(2E-5));
    end
    DiamList(n_iter) = Diam;
    OASPLList2(n_iter) = OASPL;
    n_iter = n_iter + 1;
    
end

%%


plot(RList,OASPLList,'-')

%%
plot(DiamList,OASPLList2,'-')

clc;
clear;
p = Parameters();

for n=1:p.sections
    for omega=0:100:15000
        for psi=0:pi/8:2*pi
            
        end
    end
end


% %% TRAILING EDGE
% tic;
%S_pp= @(omega) S_pp_TE_PS+S_pp_TE_SS+S_pp_LE(n);
%P_rms_BB(n)=sqrt(integral(S_pp(n),-inf,inf));
%P_rms_total_BB = sum(P_rms_BB);
%P_rms_total_BB_dB = 20.*log10(P_rms_total./p.pref);

freq= @(psi,freq_nondop) freq_nondop*(1+p.omega*p.R1*sin(theta)*sin(psi)/sqrt(p.cpcv*p.Rgas*p.temp*(1-p.Mach^2*sin(theta)^2)))
%% EXTRA
% LEADING EDGE
% 
% %1.- Radiation integral
% 
% k_bar(n) = @(freq) freq.*p.chord(n)./(2.*p.c);
% 
% epsilon(n)=1./sqrt(1+1./(4.*k_bar(n)));
% sigma=sqrt(p.x.^2+(1-p.Mach.^2).*(p.y.^2+p.z.^2));
% k_bar_y(n)=p.k_y.*p.chord(n)./2; %k_y might be a vector
% k_bar_x(n)=p.k_x.*p.chord(n)./2; %k_x might be a vector
% mu_bar(n)=@(freq) freq.*p.chord(n)./(2.*p.c.*p.beta.^2); %beta might be a vector 
% B(n)=k_bar_x(n)-p.Mach.*mu_bar(n)+k_bar(n);
% C(n)=k_bar_x(n)-mu_bar(n).*(p.x./sigma-p.Mach);
% theta(n)=k_bar(n)-mu_bar(n).*p.x./sigma;
% 
% 
% 
% %2.- Wall-pressure wave-number frequency spectrum
% b = 1./2.1; %Calibration constant for l_y
% kappa=0.41;
% 
% %SS
% alpha_SS(n)=p.V./p.U_c; %p.U_c might be a vector (for sure)
% H_SS(n)=((1+1i).*exp(-4.*1i.*k_bar(n)).*(1-theta(n).^2))./(2.*sqrt(pi).*...
%     (alpha_SS(n)-1).*p.k_x.*sqrt(B(n)));%k_x might be a vector
% G_SS(n)=(1+epsilon(n)).*exp(1i.*(2.*k_bar(n)+theta(n))).*sin(theta(n)-2.*k_bar(n))./...
%     (theta(n)-2.*k_bar(n))+(1-epsilon(n)).*exp(1i.*(-2.*k_bar(n)+theta(n))).*...
%     sin(theta(n)+2.*k_bar(n))./(theta(n)+2.*k_bar(n))+(1+epsilon(n)).*(1-1i)./...
%     (2.*(theta(n)-2.*k_bar(n))).*exp(4.*1i.*k_bar(n)).*fresnel(4.*k_bar(n))-(1-epsilon(n))...
%     .*(1+1i)./(2.*(theta(n)+2.*k_bar(n))).*exp(4.*1i.*k_bar(n)).*fresnel(4.*k_bar(n))...
%     +exp(2.*1i.*theta(n))./2.*sqrt(2.*k_bar(n)./theta(n)).*fresnel(2.*theta(n)).*...
%     ((1-epsilon(n)).*(1+1i)./(theta(n)+2.*k_bar(n))-(1+epsilon(n)).*(1+1i)./(theta(n)-2.*k_bar(n))); %I use Fresnel cosine func
% 
% I_TE_1_SS(n)=1i.*exp(2.*1i.*C(n))./C(n).*((1+1i).*exp(-2.*1i.*C(n)).*sqrt(B(n)./(B(n)-C(n))).*fresnel(2.*B(n)-2.*C(n))-(1+1i).*fresnel(2.*B(n))+1); %I use the Fresnel integral of cosine, ask if we should use the sine one.
% 
% Epsilon_param(n)=(exp(4.*1i.*k_bar(n)).*(1-(1+1i).*fresnel(4.*k_bar(n))))
% 
% Corr_epsilon(n)=Epsilon_param(n)+(epsilon(n)-1)*imag(Epsilon_param(n));
% 
% I_TE_2_SS(n)=H_SS(n).*Corr_epsilon(n)+H_SS(n).*(exp(-2.*1i.*theta(n))+1i.*(theta(n)+k_bar_x(n)+p.Mach.*mu_bar(n)-k_bar(n)).*G_SS(n)); %Take a look at definition of epsilon (slide 8 lecture 4), it should not be an exponent.
% 
% I_TE_SS(n)=I_TE_1_SS(n)+I_TE_2_SS(n);
% 
% Re_T_SS(n)=p.u_tau.*p.delta_SS(n).*sqrt(p.C_f./2)./p.kinvis; %p.u_tau, p.C_f./2 might be array
% beta_c_SS(n)=theta(n)./p.tau_w_SS(n).*p.dP_dx_SS(n);
% fun_PI_SS=@(PI_SS_var) 2.*PI_SS_var-log(1+PI_SS_var)-kappa.*p.V./...
%     p.u_tau-log(p.deltastar_SS(n).*p.V./p.kinvis)-0.51.*kappa-log(kappa);
% PI_SS(n) = 2;%fzero(fun_PI_SS,1);
% 
% phi_pp_omega_SS(n)=@(freq) (p.tau_w_SS(n).^2.*p.deltastar_SS(n)./(p.V).*0.78.*...
%     (1.8.*PI_SS(n).*beta_c_SS(n)+6).*(freq*p.deltastar_SS(n)./p.V).^2)./...
%     (((freq.*p.deltastar_SS(n)./p.V).^0.75+0.105).^3.7+...
%     (3.76.*Re_T_SS(n).^(-0.57).*(freq.*p.deltastar_SS(n)./p.V)).^7);
% 
% l_y_TE_SS(n)=@(freq) b.*p.U_c./freq; %p.U_c
% phi_pp_SS(n)=phi_pp_omega_SS(n).*l_y_TE_SS(n)./PI_SS(n);
% 
% intS_PP_TE_SS= @(k_bar_y) (phi_pp_SS(n).*sin(p.R1./p.chord(n).*(k_bar_y-k_bar(n).*p.z./sigma)).^2./...
%     ((k_bar_y-k_bar(n).*p.z./sigma).^2).*abs(I_TE_SS(n)).^2);
% 
% S_pp_TE_SS(n)=(k_bar(n).*p.z./(2.*PI_SS(n).*sigma.^2)).^2.*2.*p.chord(n).*...
%     integral(intS_PP_TE_SS,-inf,inf,'ArrayValued',1);
% 
% %PS
% alpha_PS(n)=p.V./p.U_c; %p.U_c might be a vector (for sure)
% H_PS(n)=((1+1i).*exp(-4.*1i.*k_bar(n)).*(1-theta(n).^2))./(2.*sqrt(pi).*...
%     (alpha_PS(n)-1).*p.k_x.*sqrt(B(n)));%k_x might be a vector
% G_PS(n)=(1+epsilon(n)).*exp(1i.*(2.*k_bar(n)+theta(n))).*sin(theta(n)-2.*k_bar(n))./...
%     (theta(n)-2.*k_bar(n))+(1-epsilon(n)).*exp(1i.*(-2.*k_bar(n)+theta(n))).*...
%     sin(theta(n)+2.*k_bar(n))./(theta(n)+2.*k_bar(n))+(1+epsilon(n)).*(1-1i)./...
%     (2.*(theta(n)-2.*k_bar(n))).*exp(4.*1i.*k_bar(n)).*fresnel(4.*k_bar(n))-(1-epsilon(n))...
%     .*(1+1i)./(2.*(theta(n)+2.*k_bar(n))).*exp(4.*1i.*k_bar(n)).*fresnel(4.*k_bar(n))...
%     +exp(2.*1i.*theta(n))./2.*sqrt(2.*k_bar(n)./theta(n)).*fresnel(2.*theta(n)).*...
%     ((1-epsilon(n)).*(1+1i)./(theta(n)+2.*k_bar(n))-(1+epsilon(n)).*(1+1i)./(theta(n)-2.*k_bar(n))); %I use Fresnel cosine func
% 
% I_TE_1_PS(n)=1i.*exp(2.*1i.*C(n))./C(n).*((1+1i).*exp(-2.*1i.*C(n)).*sqrt(B(n)./(B(n)-C(n))).*fresnel(2.*B(n)-2.*C(n))-(1+1i).*fresnel(2.*B(n))+1); %I use the Fresnel integral of cosine, ask if we should use the sine one.
% I_TE_2_PS(n)=H_PS(n).*(exp(4.*1i.*k_bar(n)).*(1-(1+1i).*fresnel(4.*k_bar(n)))).^epsilon(n)+H_PS(n).*(exp(-2.*1i.*theta(n))+1i.*(theta(n)+k_bar_x(n)+p.Mach.*mu_bar(n)-k_bar(n)).*G_PS(n)); %Take a look at definition of epsilon (slide 8 lecture 4), it should not be an exponent.
% 
% I_TE_PS(n)=I_TE_1_PS(n)+I_TE_2_PS(n);
% 
% Re_T_PS(n)=p.u_tau.*p.delta_PS(n).*sqrt(p.C_f./2)./p.kinvis; %p.u_tau, p.C_f./2 might be array
% beta_c_PS(n)=theta(n)./p.tau_w_PS(n).*p.dP_dx_PS(n);
% fun_PI_PS = @(PI_PS_var) (2.*PI_PS_var-log(1+PI_PS_var)-kappa.*p.V./...
%     p.u_tau-log(p.deltastar_PS(n).*p.V./p.kinvis)-0.51.*kappa-log(kappa));
% PI_PS(n)= 2; %fzero(fun_PI_PS,1);
% phi_pp_omega_PS(n)=@(freq) (p.tau_w_PS(n).^2.*p.deltastar_PS(n)./(p.V).*0.78.*...
%     (1.8.*PI_PS(n).*beta_c_PS(n)+6).*(freq.*p.deltastar_PS(n)./p.V).^2)./...
%     (((freq.*p.deltastar_PS(n)./p.V).^0.75+0.105).^3.7+(3.76.*...
%     Re_T_PS(n).^(-0.57).*(freq.*p.deltastar_PS(n)./p.V)).^7);
% 
% l_y_TE_PS(n)=@(freq) b.*p.U_c./freq; %p.U_c
% phi_pp_PS(n)= phi_pp_omega_PS(n).*l_y_TE_PS(n)./PI_PS(n);
% intS_pp_TE_PS = @(k_bar_y) phi_pp_PS(n).*sin(p.R1./p.chord(n).*...
%     (k_bar_y-k_bar(n).*p.z./sigma)).^2./((k_bar_y-k_bar(n).*p.z./sigma).^2)...
%     .*abs(I_TE_PS(n)).^2;
% S_pp_TE_PS(n)=(k_bar(n).*p.z./(2.*PI_PS(n).*sigma.^2)).^2.*2.*...
%     p.chord(n).*integral(intS_pp_TE_PS,-inf,inf,'ArrayValued',1);
% 
% 
% S_pp_TE(n) = S_pp_TE_SS(n) + S_pp_TE_PS(n);
% toc;
% end
% %S_pp= @(omega) S_pp_TE_PS+S_pp_TE_SS+S_pp_LE(n);
% %P_rms_BB(n)=sqrt(integral(S_pp(n),-inf,inf));
% %P_rms_total_BB = sum(P_rms_BB);
% %P_rms_total_BB_dB = 20.*log10(P_rms_total./p.pref);
% 
% freq= @(psi,freq_nondop) freq_nondop*(1+p.omega*p.R1*sin(theta)*sin(psi)/sqrt(p.cpcv*p.Rgas*p.temp*(1-p.Mach^2*sin(theta)^2))
% %% EXTRA
% % LEADING EDGE
% % 
% % %1.- Radiation integral
% % theta_plus(n)=k_bar(n)+mu_bar(n)./sigma;
% % theta_minus(n)=k_bar(n)-mu_bar(n)./sigma;
% % theta_bar(n)=k_bar_x(n).*(1-p.Mach.*p.x./sigma)./(1-p.Mach.^2)-pi./4;
% % I_LE_1(n)=-1./pi.*sqrt(2./(theta_minus(n).*(k_bar_x(n)+(1-p.Mach.^2).*...
% %     k_bar(n).*theta_minus(n)))).*exp(-1i.*theta(n)).*fresnel(2.*theta(n));
% % I_LE_2(n)=exp(-1i.*theta_bar(n))./(pi.*theta(n).*sqrt(2.*pi.*(k_bar_x(n)+...
% %     (1-p.Mach.^2).*k_bar(n)))).*(1i.*(1-exp(2.*1i.*theta_minus(n)))-...
% %     (1+1i).*(fresnel(4.*k_bar(n))-exp(2.*1i.*theta_minus(n)).*...
% %     sqrt(2.*k_bar(n)./theta_plus(n)).*fresnel(2.*theta_plus(n))));
% % 
% % I_LE(n)=I_LE_1(n)+I_LE_2(n);
% % 
% % %2.- Phi omega omega
% % L_t=0.4.*p.turb_kin_en.^1.5./p.diss_rate; %check if they change (maybe not)
% % u_bar_prime=sqrt(2.*p.turb_kin_en./3);
% % k_e=sqrt(pi)./L_t.*gamma(5./6)./gamma(1./3);
% % phi_omegaomega_omega = u_bar_prime.^2.*L_t./(2.*pi.*p.V).*...
% %     (1+8./3.*(p.k_x+k_e).^2)./((1+(p.k_x./k_e).^2).^(211./6)); %check p.k_x
% % l_y_LE=(gamma(1./3)./gamma(5./6)).^2.*(p.k_x./k_e).^2./...
% %     ((3+8.*(p.k_x+k_e).^2).*sqrt(1+(p.k_x+k_e).^2));
% % 
% % phi_omegaomega=p.V.*phi_omegaomega_omega.*l_y_LE./pi;
% % 
% % %3.- Power spectral density
% % intS_pp_LE = @(k_y_LE) (phi_omegaomega.*sin(p.R1./2.*(k_y_LE./sigma-k_y_LE)).^2./...
% %     (PI_SS(n)./2.*(k_y_LE./sigma-k_y_LE).^2.*abs(I_LE(n)).^2));
% % 
% % S_pp_LE=(p.density.*k_bar(n).*p.z./(sigma.^2)).^2.*pi.*p.V.*p.R1./2.*... %check if k_y array
% %     integral(intS_pp_LE,-inf,inf);
% % 
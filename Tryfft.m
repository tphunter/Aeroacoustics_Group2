t = linspace(0,0.1,1000);
%0:1/50:10-1/50;                     
x = sin(2*pi*100*t);% + sin(2*pi*20*t)+1;
% plot(t,x)
y = fft(x);     
f = (0:length(y)-1)*200/length(y);
% plot(f,abs(y))
title('Magnitude')
n = length(x);                         
fshift = (-n/2:n/2-1)*(250/n);
yshift = fftshift(y);
% plot(fshift,abs(yshift))


Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:L-1)*T;        % Time vector
S = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t);
plot(1000*t,S)
title('Signal')
xlabel('t (milliseconds)')
ylabel('X(t)')

X = S + 2*randn(size(t));
plot(1000*t(1:50),X(1:50))
title('Signal Corrupted with Zero-Mean Random Noise')
xlabel('t (milliseconds)')
ylabel('X(t)')
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

Y = fft(S);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L-1))/L;
plot(f,P2) 
title('Single-Sided Amplitude Spectrum of S(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')




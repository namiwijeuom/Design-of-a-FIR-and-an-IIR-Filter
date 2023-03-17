% Index number - 200733D
% According to the index number,
clc;

close all;
A=7;
B=3;
C=3; 

% Design Specifications
Ap  =   0.1+(0.01* A);          % Maximum passband ripple 
Aa  =   50 + B;                 % Minimum stopband attenuation 
wp1 =   (C * 100) + 400;        % Lower passband edge 
wp2 =   (C * 100) + 900;        % Upper passband edge 
wa1 =   (C * 100) + 100;        % Lower stopband edge 
wa2 =   (C * 100) + 1100;        % Upper stopband edge 
ws  =    2*(( C * 100) + 1500) ; % Sampling frequency

% Sampling period
Ts = 2*pi/ws;

% Sampling frequency
fs = 1/Ts;

% Passband and stopband frequencies
wp = [wp1 wp2];
ws1 = [wa1 wa2];

% Apply pre - warping transformation
wp_warped = 2*tan(wp*Ts/2)/Ts;
ws_warped = 2*tan(ws1*Ts/2)/Ts;

%Getting the order and the cutoff frequencies 
[n,Wp] = ellipord(wp_warped, ws_warped , Ap, Aa,"s");  
disp(n);
disp(Wp);


[num,den] = ellip(n,Ap,Aa,Wp,"bandpass","s");

%Applying Bilinear transform
[dnum,dden]=bilinear(num,den,fs);

%Getting the required digital filter coefficients
disp(dnum);
disp(dden);
printsys (dnum , dden )

[H,w]=freqz(dnum,dden,2001);

Hdb = 20*log10(abs(H));

%Plotting the Magnitude response of the filter
figure;
plot([flip(-w); w], [flip(Hdb); Hdb])
xlabel('\Omega (rad/sample)')
ylabel('Magnitude (dB)')
title('Magnitude response')
ax = gca;
ax.YLim = [-200 20];
ax.XLim = [-pi pi];
grid on;
grid minor;

%Plotting the Magnitude response in the Passband
figure;
plot(w, Hdb);
xlabel('\Omega (rad/sample)')
ylabel('Magnitude (dB)')
title('Magnitude response in passband')
ax = gca;
ax.YLim = [-0.2 0.2];
ax.XLim = [wp1*Ts wp2*Ts];
grid on;
grid minor;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Testing
% 
% %Input Signal Generation
% no_of_samples = 300;
% w1 = wa1/2;             %Middle freueny of the lower stopband
% w2 = (wp1 + wp2)/2;     %Middle frequency of the passband
% w3 = (wa2 + ws/2)/2;    %Middle frequency of the upper stopband
% 
% Input_signal = zeros(no_of_samples,1);
% 
% for i = 1:no_of_samples
%     Input_signal(i) = sin(w1*i*Ts) + sin(w2*i*Ts) + sin(w3*i*Ts);
% end
% 
% %Plotting the Input Signal
% figure;
% subplot(2,1,1);
% stem(Input_signal);
% title('Input Signal');
% xlabel('Samples');
% ylabel('Amplitude');
% 
% %Actual Output Signal
% lenfft = length(Input_signal) + length(H) - 1;
% 
% XF = fft(Input_signal,lenfft);  %Taking the FT of the input signal
% HF = fft(H,lenfft);       %Taking the FT of the filter
% YF = HF.*XF;                
% 
% Actual_Output = ifft(YF,lenfft); %Taking the Inverse FT of YF
% 
% subplot(2,1,2);
% plot(Actual_Output);
% axis([0 no_of_samples -2 2]);
% title('Actual Output Signal');
% xlabel('Samples');
% ylabel('Amplitude');

% Index number - 200733D

% According to the index number,
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

%Step 1 - Selecting the transistion bandwidth and cutoff points
Bt = min(wp1 - wa1, wa2 - wp2); %Transistion bandwidth is the minumum of lower transistion width and upper transistion width
wc1 = wp1 - Bt/2 ; %Lower cutoff frequency
wc2 = wp2 + Bt/2 ; %Upper cutoff frequency
T = 2*pi/ws; %Sampling period
disp(Bt);disp(wc1);disp(wc2);disp(T);

%Step 2 - Choosing Delta
delta_p = (10^( Ap/20 ) - 1) / (10^( Ap/20 ) + 1);
delta_a = 10^( -Aa/20 );
delta = min(delta_p,delta_a);
disp(delta)

%Step 3 - Finding the actual stopband attenuation for the defined delta
actual_Aa = -20*log10(delta);

%Step 4 - Choosing the parameter Alpha based on the actual stopband attenuation
if (actual_Aa <= 21)
    alpha = 0;
elseif (actual_Aa <= 50)
    alpha = 0.5842*(actual_Aa - 21)^0.4 + 0.07886*(actual_Aa - 21);
else
    alpha = 0.1102*(actual_Aa - 8.7);
end
disp(alpha)

%Step 5 - Choosing the parameter D based on the actual stop band attenuation
if (actual_Aa <= 21)
    D = 0.9222;
else
    D = (actual_Aa - 7.95) / 14.36;
end
disp(D)

%Selecting the lowest off value N which satifies the inequality N >= ( (ws*D) / Bt )+1
N = ceil(( ws*D / Bt )+1);
if(mod(N,2)==0)
    N = N + 1;
end

disp("The order of the filter is ")
disp(N)

%Forming Wk(nT), the Kaiser Window function 
wk = zeros(N,1);
for n = -(N - 1)/2:(N - 1)/2
   beta = alpha * (1 - (2*n/(N - 1))^2)^0.5;
   numerator = my_Bessel_func(beta);
   denominator = my_Bessel_func(alpha);
   wk(n+(N-1)/2+1) = numerator/denominator;
end

stem(wk,'filled');
title('Kaiser Window Function');
xlabel('n');
ylabel('Wn');
grid on;

%Step 6 - Getting h(nT) using the formed wk(Nt)
h = zeros(N,1);
h(38) = (2/ws)*(wc2 - wc1);
for n = -(N-1)/2:(N-1)/2
   if n==0
      h(n+(N-1)/2+1) = (2/ws)*(wc2 - wc1);
   else
      h(n+(N-1)/2+1) = (1/(n*pi)) * (sin(wc2*n*T) - sin(wc1*n*T));
   end
end
figure;
stem(h,'filled');
title('h(nT)');

digital_filter = h.*wk;

%Impulse response of the Filter
figure;
stem(digital_filter,'filled');
title('Impulse Response of the Digital Filter');
xlabel('Samples (which are given by n + (N-1)/2)');
ylabel('Amplitude');
grid on;

% Magnitude Response of the Filter in rad/s
[h,w]=freqz(digital_filter);
magnitude = 20*log10(abs(h));
analogFreq = w*ws/(2*pi);
figure;
plot(analogFreq,magnitude);
grid ON;
title('Magnitude Response of the Digital Filter');
xlabel('Angular Frequency in rad/s');
ylabel('Magnitude(dB)');

% Magnitude Response of the Filter in rad/sample
fvtool(digital_filter);

%Magnitude Response in Passband
[amp, digiFreq] = freqz(digital_filter);
analogFreq = digiFreq*ws/(2*pi);
ampdb = 20*log10(abs(amp));
figure;
plot(analogFreq, ampdb);
axis([wp1 wp2 -0.05 0.05]);
title('Magnitude Response in Passband');
xlabel('Angular Frequency in rad/s');
ylabel('Magnitude(dB)');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Testing

%Input Signal Generation
no_of_samples = 300;
w1 = wa1/2;             %Middle freueny of the lower stopband
w2 = (wp1 + wp2)/2;     %Middle frequency of the passband
w3 = (wa2 + ws/2)/2;    %Middle frequency of the upper stopband

Input_signal = zeros(no_of_samples,1);

for i = 1:no_of_samples
    Input_signal(i) = sin(w1*i*T) + sin(w2*i*T) + sin(w3*i*T);
end

%Plotting the Input Signal
figure;
subplot(3,1,1);
stem(Input_signal);
title('Input Signal');
xlabel('Samples');
ylabel('Amplitude');

%Expected Output Signal
expected_output = zeros(no_of_samples,1);

for i = 1:no_of_samples
    expected_output(i) = sin(w2*i*T);
end

subplot(3,1,2);
stem(expected_output);
axis([0 no_of_samples -2 2]);
title('Expected Output Signal');
xlabel('Samples');
ylabel('Amplitude');

%Actual Output Signal
lenfft = length(Input_signal) + length(digital_filter) - 1;

XF = fft(Input_signal,lenfft);  %Taking the FT of the input signal
HF = fft(digital_filter,lenfft);       %Taking the FT of the filter
YF = HF.*XF;                

Actual_Output = ifft(YF,lenfft); %Taking the Inverse FT of YF

subplot(3,1,3);
stem(Actual_Output);
axis([0 no_of_samples -2 2]);
title('Actual Output Signal');
xlabel('Samples');
ylabel('Amplitude');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Bessel function which is required to find the Kaiser window function
function [ result ] = my_Bessel_func( x )
    j = 1;
    result = 0;
    term = 10;
    while (term > 10^(-6))
        term = (((x/2)^j)/(factorial(j)))^2;
        result = result + term;
        j = j + 1;
    end
    
    result = result + 1;
end
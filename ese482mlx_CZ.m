%ese482 matlab project
%Claire Zhang 
%12/09/2018

echo on

% This MATLAB file designs two zero phase bandpass filters - one elliptic filter and the other optimal filter in the min-max sense, that can meet the secifications listed below:
% pass band(3dB down) from f1=50 to f2=100
% ripple in the passband less than r=3dB
% stopband s=25dB down
% Figures are plotted to demonstrate that the designed filters meet the requriements. 
% Then the filters are applied to the provided ECG data which are sampled at 1000Hz. 
% The implementation times of the two methods are compared to pick the filter that has fewer computations, i.e, shorter implementation time. 
% This faster filter is applied to ECG data to compare the input and output signal.
% Random noise that is p percent of the maximum ECG value is added to the data set which is then filtered to compare the output result when the original data have noise and when the data do not.
% Precision of the coefficients for the faster filter implementation is reduced until the specifications can no longer be met. 
% The lowest number of bits possible to retain the filter feature is calculated.
pause;


% elliptic filter
f1=50; f2=100; n = 4; Rp = 3; Rs=25; Wp=[f1/500 f2/500]; % 500 is half the sampling rate
[b,a] = ellip(n,Rp,Rs,Wp);
figure
freqz(b,a)
title('Frequency Response of Elliptic Filter')
% Figure shows that the specifications are satisfied.
pause;
figure
zplane(b,a)
title('Pole-Zero Plot of Elliptic Filter')
pause;


% optimal filter in min-max sense 
fx = [0 0.08 50/500 100/500 0.22 1];
fy = [0 0 1 1 0 0];
figure
plot(fx,fy,'linewidth',2)
xlabel('Normalized Frequency (\times \pi rad/sample)')
ylabel('Magnitude')
title('Comparison between Ideal and Designed Filter in Min-Max Sense')
hold on
coeff =firpm(113,fx,fy);
[h,w] = freqz(coeff);
plot(w/pi,abs(h),'linewidth',2)
legend('Ideal Filter','firpm')
%The plot compares the gain of the ideal filter and the firpm filter.
pause;
figure
freqz(coeff)
title('Firpm Filter Design Frequency Response')
%Figure shows that the firpm filter also meets the specifications.
pause;


% Compare implementation time of elliptic and firpm filters
% Elliptic filter is implemented using function filtfilt in order to perform zero-phase filtering.
% time ellipic filter implementation
echo off
load ar0922d;
tic
    for i=1:20000
            rmadip_filt = filtfilt(b,a,rmadip);
    end
timeElip = toc;
echo on
timeElip = timeElip/20000 % Implementation time of elliptic filter
pause;
% Firpm filter is implemented by computing discrete fourier transform using fast fourier transform (FFT) algorithm.
% Mathematically, instead of directly calcualting the convolution sum of the filter and input data, FFT is first applied to both the filter and input data. 
% The results are then multiplied together and converted back to the convolution sum by using inverse FFT.
% time firpm filter implementation
echo off
rmadip_new = zeros(3,588);
rmadip_new_filt = zeros(3,588);
rmadip_new_fft = zeros(3,588);
for i=1:3
    for j=1:588
    rmadip_new(i,j) = rmadip(j,i); %change each column to a row
    end
end
tic 
    for j=1:20000
        filter_fft = fft(coeff,length(rmadip_new(1,:))); % fft of firpm filter
        for i=1:3
            rmadip_new_fft(i,:)=fft(rmadip_new(i,:));
            rmadip_new_filt(i,:)=ifft(filter_fft.*rmadip_new_fft(i,:));
        end
    end
timefirpm = toc;
echo on
timefirpm = timefirpm/20000 % implementation time of firpm filter
pause;
filtertype = {'Elliptic Filter';'FIRPM'};
implementation_time = [timeElip;timefirpm];
T = table(filtertype,implementation_time)
% The firpm filter is a faster implementation. 
% The elliptic filter is a slower implementation because it uses function filtfilt to avoid phase distortion. 
% The filtfilt function applies both forward and backward filtering to the data, resulting in an implementation time that is about doubled because it calculates the convolution sum twice.
% The firpm itself is in linear phase. But the phase is corrected by a direct time delay shift, which is a faster mathematical operation. 
% Also since the implementation of firpm does not go through convolution sum, the calculation itself is easier.
% Since firpm is the faster implementation, it is used for the rest of this MATLAB project.
pause;


% compare input and output ECG signals
% ECG data from each of the three electrodes are plotted separately in three subplots.
echo off
N = length(coeff);
delay = N/2+1; 
% find filtered data using firpm and correct time delay by using a direct shift
for i=1:3
    rmadip_new_filt_shifted(i,:)= rmadip_new_filt(i,delay+1:end);
end
figure
subplot(3,1,1);plot(rmadip_new(1,:),'linewidth',2);hold on
plot(rmadip_new_filt_shifted(1,:),'linewidth',2);hold off
xlabel('Sample')
ylabel('Magnitude')
title('Comparison between Input and Shifted FIRPM Output of X Lead ECG Data')
subplot(3,1,2);plot(rmadip_new(2,:),'linewidth',2);hold on
plot(rmadip_new_filt_shifted(2,:),'linewidth',2);hold off
xlabel('Sample')
ylabel('Magnitude')
title('Comparison between Input and Shifted FIRPM Output of Y Lead ECG Data')
subplot(3,1,3);plot(rmadip_new(3,:),'linewidth',2);hold on
plot(rmadip_new_filt_shifted(3,:),'linewidth',2);hold off
xlabel('Sample')
ylabel('Magnitude')
title('Comparison between Input and Shifted FIRPM Output of Z Lead ECG Data')
echo on
% The filtered data have lowered amplitude.
% The frequency range in which the abnormity in ECG data is present stands out in the filtered data.
pause;


% Add noise to ECG data and compare input and output
% p value is found through trials and errors
echo off
p = 0.0062;
for i=1:3
    B = randn(1,length(rmadip_new(1,:)));
    B_normalized = (B./max(B)).*(p*max(rmadip_new(i,:)));
    rmadip_new_noise(i,:) = rmadip_new(i,:)+B_normalized;
end
 for i=1:3
    rmadip_new_noise_fft(i,:)=fft(rmadip_new_noise(i,:));
    rmadip_new_noise_filt(i,:)=ifft(filter_fft.*rmadip_new_noise_fft(i,:));
 end
% compare input and outout of the signal with noise 
for i=1:3
    rmadip_new_noise_filt_shifted(i,:)= rmadip_new_noise_filt(i,delay+1:end);
end
figure
subplot(3,1,1);plot(rmadip_new_noise_filt_shifted(1,:),'linewidth',2);hold on
plot(rmadip_new_filt_shifted(1,:),'linewidth',2);hold off
xlabel('Sample')
ylabel('Magnitude')
title('Comparison between Shifted FIRPM Output of X Lead ECG Data with and without Noise')
legend('With Noise','Without Noise');
subplot(3,1,2);plot(rmadip_new_noise_filt_shifted(2,:),'linewidth',2);hold on
plot(rmadip_new_filt_shifted(2,:),'linewidth',2);hold off
xlabel('Sample')
ylabel('Magnitude')
title('Comparison between Shifted FIRPM Output of Y Lead ECG Data with and without Noise')
legend('With Noise','Without Noise');
subplot(3,1,3);plot(rmadip_new_noise_filt_shifted(3,:),'linewidth',2);hold on
plot(rmadip_new_filt_shifted(3,:),'linewidth',2);hold off
xlabel('Sample')
ylabel('Magnitude')
legend('With Noise','Without Noise');
title('Comparison between Shifted FIRPM Output of Z Lead ECG Data with and without Noise')
xlim([0,530])
echo on
% The largest p value at which noise is not observed in the filtered result is p=0.62%.
pause;


% calculate rms error 
% rms error is calculated based on the following mathematical equations:
% rms error = sqrt((1/n)*(sum(difference each entry in the two data sets)^2))
echo off
difference =rmadip_new_filt_shifted-rmadip_new_noise_filt_shifted;
rmsvalue = 0;
for i=1:3
    for j=1:length(difference(1,:))
        rmsvalue = rmsvalue + (difference(i,j))^2;
    end
end
echo on
rmsvalue = sqrt(rmsvalue/(3*length(difference(1,:))))
% The rms error is found to be 0.0511.


% calculate the bits 
% number of bits is found through trials and errors
n = 10;
coeff_new = round(coeff.*(2^n))./2^n;
freqz(coeff_new)
title('Frequency Response of the Reduced FIRPM Filter')
% When the number of bits goes down to about 10, the firpm filter features are still preserved and the specifications satisfied. 
% But the number of bits cannot go under 10, below which the filter is getting unable to meet the specifications. 
% For example, when the number of bits goes below 10, the stop band can not long be consistently 25dB down
% The lowest possible bits is found to be 10.
echo off



 


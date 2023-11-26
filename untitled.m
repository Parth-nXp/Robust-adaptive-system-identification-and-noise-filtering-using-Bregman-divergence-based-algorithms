clc;
clear all;
close all;

signal = audioread('handel.wav');
signal = signal + 3*(abs(min(signal)))*ones(size(signal));
figure;
subplot(5,2,1)
plot(signal)
title("Original Signal")
mu_KLLMS = 0.02;
mu_ISLMS = 0.02;
mu_AALMS = 0.02;
mu_BLMS = 0.02;
channel_taps = 4;
experiment = 1;
alpha = 4;
beta = 2; 
RandStream.setGlobalStream(RandStream('mt19937ar','seed',109));
filter_weights = randn(channel_taps,1)+10e-6*ones(channel_taps,1);
MSE_KLLMS_main = zeros(size(signal));
MSE_ISLMS_main = zeros(size(signal));
MSE_AALMS_main = zeros(size(signal));
MSE_BLMS_main = zeros(size(signal));
filter_output_KLLMS_main = zeros(size(signal));
filter_output_ISLMS_main = zeros(size(signal));
filter_output_AALMS_main = zeros(size(signal));
filter_output_BLMS_main = zeros(size(signal));
for i = 1:experiment
    noise = sqrt(0.05)*randn(1,length(signal));
    noisy_signal = signal' + noise;
    w_KLLMS = filter_weights;
    w_ISLMS = filter_weights;
    w_AALMS = filter_weights;
    w_BLMS = filter_weights;
    filter_output_KLLMS = zeros(size(signal));
    filter_output_ISLMS = zeros(size(signal));
    filter_output_AALMS = zeros(size(signal));
    filter_output_BLMS = zeros(size(signal));
    error_KLLMS = zeros(size(signal));
    error_ISLMS = zeros(size(signal));
    error_AALMS = zeros(size(signal));
    error_BLMS = zeros(size(signal));
    MSE_KLLMS = zeros(size(signal));
    MSE_ISLMS = zeros(size(signal));
    MSE_AALMS = zeros(size(signal));
    MSE_BLMS = zeros(size(signal));
    un = zeros(channel_taps,1);
    for n=1:length(signal)
        un = [noisy_signal(n); un(1:end-1)];
        filter_output_KLLMS(n) = un'*w_KLLMS;
        filter_output_ISLMS(n) = un'*w_ISLMS;
        filter_output_AALMS(n) = un'*w_AALMS;
        filter_output_BLMS(n) = un'*w_BLMS;
        error_KLLMS(n) = signal(n) - filter_output_KLLMS(n);
        error_ISLMS(n) = signal(n) - filter_output_ISLMS(n);
        error_AALMS(n) = signal(n) - filter_output_AALMS(n);
        error_BLMS(n) = signal(n) - filter_output_BLMS(n);
        w_KLLMS = w_KLLMS + mu_KLLMS * un *((signal(n)/filter_output_KLLMS(n))-1);
        w_ISLMS = w_ISLMS + mu_ISLMS * un * ((signal(n)/(filter_output_ISLMS(n))^2)-(1/(filter_output_ISLMS(n))));
        w_AALMS = w_AALMS + (mu_AALMS/alpha) * un * ((signal(n)^alpha/(filter_output_AALMS(n))^alpha)-1);
        w_BLMS = w_BLMS + mu_BLMS * un * (filter_output_BLMS(n))^(beta-1)*error_BLMS(n);
        MSE_KLLMS(n) = error_KLLMS(n)^2;
        MSE_ISLMS(n) = error_ISLMS(n)^2;
        MSE_AALMS(n) = error_AALMS(n)^2;
        MSE_BLMS(n) = error_BLMS(n)^2;
    end
    filter_output_KLLMS_main = filter_output_KLLMS_main + filter_output_KLLMS;
    filter_output_ISLMS_main = filter_output_ISLMS_main + filter_output_ISLMS;
    filter_output_AALMS_main = filter_output_AALMS_main + filter_output_AALMS;
    filter_output_BLMS_main = filter_output_BLMS_main + filter_output_BLMS;
    MSE_KLLMS_main = MSE_KLLMS_main + MSE_KLLMS;
    MSE_ISLMS_main = MSE_ISLMS_main + MSE_ISLMS;
    MSE_AALMS_main = MSE_AALMS_main + MSE_AALMS;
    MSE_BLMS_main = MSE_BLMS_main + MSE_BLMS;
end
subplot(5,2,2)
plot(noisy_signal)
title("Noisy Signal")
MSE_KLLMS_main = MSE_KLLMS_main/experiment;
MSE_ISLMS_main = MSE_ISLMS_main/experiment;
MSE_AALMS_main = MSE_AALMS_main/experiment;
MSE_BLMS_main = MSE_BLMS_main/experiment;
MSE_KLLMS_main = MSE_KLLMS_main/max(MSE_KLLMS_main);
MSE_ISLMS_main = MSE_ISLMS_main/max(MSE_ISLMS_main);
MSE_AALMS_main = MSE_AALMS_main/max(MSE_AALMS_main);
MSE_BLMS_main = MSE_BLMS_main/max(MSE_BLMS_main);
filter_output_KLLMS_main = filter_output_KLLMS_main/experiment;
filter_output_ISLMS_main = filter_output_ISLMS_main/experiment;
filter_output_AALMS_main = filter_output_AALMS_main/experiment;
filter_output_BLMS_main = filter_output_BLMS_main/experiment;
subplot(5,2,3)
plot(filter_output_KLLMS_main)
% title("KLLMS recovered Signal")
ylim([1.5 3.5])
subplot(5,2,4)
plot(10*log10(MSE_KLLMS_main));
title("KLLMS MSE")
subplot(5,2,5)
plot(filter_output_ISLMS_main)
% title("ISLMS recovered Signal")
ylim([1.5 3.5])
subplot(5,2,6)
plot(10*log10(MSE_ISLMS_main));
title("ISLMS MSE")
subplot(5,2,7)
plot(filter_output_AALMS_main)
% title("AALMS recovered Signal")
ylim([1.5 3.5])
subplot(5,2,8)
plot(10*log10(MSE_AALMS_main));
title("AALMS MSE")
subplot(5,2,9)
plot(filter_output_BLMS_main)
% title("BLMS recovered Signal")
ylim([1.5 3.5])
subplot(5,2,10)
plot(10*log10(MSE_BLMS_main));
title("BLMS MSE")
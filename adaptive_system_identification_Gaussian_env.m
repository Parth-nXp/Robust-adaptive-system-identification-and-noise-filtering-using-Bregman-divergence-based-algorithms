clc;
clear all;
close all;

mu_KLLMS = 0.002;
mu_ISLMS = 0.002;
mu_AALMS = 0.002;
mu_BLMS = 0.002;
mu_LMS = 0.002;
channel_taps = 6;
experiment = 1000;
alpha = 2;
beta = -2; 
iteration = 10000;
RandStream.setGlobalStream(RandStream('mt19937ar','seed',62));  %data1
filter_weights = abs(randn(channel_taps,1));
filter_weights = filter_weights/norm(filter_weights);
RandStream.setGlobalStream(RandStream('mt19937ar','seed',65));  %data1
initial_weights = abs(randn(channel_taps,1));
MSD_KLLMS_main = zeros(iteration,1);
MSD_ISLMS_main = zeros(iteration,1);
MSD_AALMS_main = zeros(iteration,1);
MSD_BLMS_main = zeros(iteration,1);
MSD_LMS_main = zeros(iteration,1);

for i = 1:experiment
    noise = sqrt(1)*randn(iteration,1);
    w_KLLMS = initial_weights;
    w_ISLMS = initial_weights;
    w_AALMS = initial_weights;
    w_BLMS = initial_weights;
    w_LMS = initial_weights;
    filter_output_KLLMS = zeros(iteration,1);
    filter_output_ISLMS = zeros(iteration,1);
    filter_output_AALMS = zeros(iteration,1);
    filter_output_BLMS = zeros(iteration,1);
    filter_output_LMS = zeros(iteration,1);
    error_KLLMS = zeros(iteration,1);
    error_ISLMS = zeros(iteration,1);
    error_AALMS = zeros(iteration,1);
    error_BLMS = zeros(iteration,1);
    error_LMS = zeros(iteration,1);
    MSD_KLLMS = zeros(iteration,1);
    MSD_ISLMS = zeros(iteration,1);
    MSD_AALMS = zeros(iteration,1);
    MSD_BLMS = zeros(iteration,1);
    MSD_LMS = zeros(iteration,1);
    un = zeros(channel_taps,1);
    for n=1:iteration
        un = [abs(normrnd(0,1)); un(1:end-1)];
        dn = un'*filter_weights + noise(n);
        filter_output_KLLMS(n) = un'*w_KLLMS;
        filter_output_ISLMS(n) = un'*w_ISLMS;
        filter_output_AALMS(n) = un'*w_AALMS;
        filter_output_BLMS(n) = un'*w_BLMS;
        filter_output_LMS(n) = un'*w_LMS;
        error_KLLMS(n) = dn - filter_output_KLLMS(n);
        error_ISLMS(n) = dn - filter_output_ISLMS(n);
        error_AALMS(n) = dn - filter_output_AALMS(n);
        error_BLMS(n) = dn - filter_output_BLMS(n);
        error_LMS(n) = dn - filter_output_LMS(n);
        w_KLLMS = w_KLLMS + mu_KLLMS * un *(dn/(filter_output_KLLMS(n))-1);
        w_ISLMS = w_ISLMS + mu_ISLMS * un * ((dn/((filter_output_ISLMS(n)+10e-6))^2)-(1/(filter_output_ISLMS(n)+10e-6)));
        w_AALMS = w_AALMS + (mu_AALMS/alpha) * un * ((dn^alpha/(filter_output_AALMS(n))^alpha)-1);
        w_BLMS = w_BLMS + mu_BLMS * un * (filter_output_BLMS(n))^(beta-1)*error_BLMS(n);
        w_LMS = w_LMS + mu_LMS * un * error_LMS(n);
        MSD_KLLMS(n) = norm(w_KLLMS-filter_weights)^2;
        MSD_ISLMS(n) = norm(w_ISLMS-filter_weights)^2;
        MSD_AALMS(n) = norm(w_AALMS-filter_weights)^2;
        MSD_BLMS(n) = norm(w_BLMS-filter_weights)^2;
        MSD_LMS(n) = norm(w_LMS - filter_weights)^2;
    end
    MSD_KLLMS_main = MSD_KLLMS_main + MSD_KLLMS;
    MSD_ISLMS_main = MSD_ISLMS_main + MSD_ISLMS;
    MSD_AALMS_main = MSD_AALMS_main + MSD_AALMS;
    MSD_BLMS_main = MSD_BLMS_main + MSD_BLMS;
    MSD_LMS_main = MSD_LMS_main + MSD_LMS;
end

MSD_KLLMS_main = MSD_KLLMS_main/experiment;
MSD_ISLMS_main = MSD_ISLMS_main/experiment;
MSD_AALMS_main = MSD_AALMS_main/experiment;
MSD_BLMS_main = MSD_BLMS_main/experiment;
MSD_LMS_main = MSD_LMS_main/experiment;

MSD_KLLMS_main = MSD_KLLMS_main/max(MSD_KLLMS_main);
MSD_ISLMS_main = MSD_ISLMS_main/max(MSD_ISLMS_main);
MSD_AALMS_main = MSD_AALMS_main/max(MSD_AALMS_main);
MSD_BLMS_main = MSD_BLMS_main/max(MSD_BLMS_main);
MSD_LMS_main = MSD_LMS_main/max(MSD_LMS_main);
figure;
plot(10*log10(MSD_KLLMS_main))
hold on
plot(10*log10(MSD_ISLMS_main))
plot(10*log10(MSD_AALMS_main))
plot(10*log10(MSD_BLMS_main))
plot(10*log10(MSD_LMS_main))
xlabel('Iteration')
ylabel('Mean Square Deviation (dB)');
legend('KLLMS','ISLMS','AALMS','BLMS','LMS')

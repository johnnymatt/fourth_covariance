%% Ex1: Define a test signal (x)
clear all;clc;close all;
% User inputs:
% f0: arbitrary constant frequency (Hz)
% Fs: sampling rate (Hz)
% L: Length of the signal
% I: Number of segments for computing FC
% snr: signal to noise ratio (added white noise) (dB)
snr = 30;
f0 = 10;
Fs = 100;
L = 100000;
I = 100;
fvec_seg = (0:floor(L / I)/2-1)*(Fs/floor(L / I)); f_spac = abs(fvec_seg(1)-fvec_seg(2));
disp(['Frequency spacing for each segment = ', num2str(f_spac), 'Hz'])
disp(['Maximum frequency in the Nyquist range = ', num2str(max(fvec_seg)), 'Hz'])

[x, t] = test_signal(f0, Fs, L);

%% Ex2: Frequency characteristics calculation of x

% Compute the two-sided spectrum using FFT
X_dft = fft(x);
X_mag2 = abs(X_dft/L);

% Compute one-sided DFT amplitude
X_mag1 = X_mag2(1:L/2+1);
X_mag1(2:end-1) = 2*X_mag1(2:end-1); % Multiply spectra in the +ve freqs by 2 as signal is real-valued

% Power Spectral Density (PSD)
X_psd = (1/(Fs*L)) * abs(X_mag1).^2;
X_psd(2:end-1) = 2*X_psd(2:end-1); % DC and Nyquist do not occur twice
X_psd_db = pow2db(X_psd);

f = 0:Fs/L:Fs/2; % Get frequency vector

% Plots
figure (1)
% Plotting the DFT
subplot(2,1,1);
plot(f, X_mag1);
title('Discrete Fourier Transform (DFT) of Test Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Plotting the PSD
subplot(2,1,2);
plot(f, X_psd_db);
title('Power Spectral Density (PSD) of Test Signal');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
grid on;

%% Ex3&4: Save .mat files
freqs_psd = [f0, 2*f0, 3*f0, 4*f0];
PSD_mags = zeros(1,length(freqs_psd));
for i = 1:length(freqs_psd)
    PSD_mags(i) = X_psd_db(f==freqs_psd(i));
end
disp(['Magnitudes for PSD are: ', num2str(PSD_mags(:)')])
freqs_fc = [f0, f0, f0, f0; f0, 2*f0, 3*f0, 4*f0; 2*f0, 2*f0, 2*f0, 2*f0];
FC_mags = zeros(size(freqs_fc,1),1);
for i = 1:size(freqs_fc,1)
    FC_mags(i) = fc(x, I, snr, Fs, freqs_fc(i,:));
end
FC_mags = abs(FC_mags); % Magnitudes
disp(['Magnitudes for FC are: ', num2str(FC_mags(:)')])

save('PSD_magnitudes.mat','PSD_mags')
save('FC_magnitudes.mat','FC_mags')
function [FC,noisy_signal] = fc(x,I,snr,Fs,fv)
    % Add white noise to the signal with a signal to noise ratio of 30dB
    noisy_signal = awgn(x, snr, 'measured');
    L_segm = floor(length(x) / I); % Length of each segment
    dft_segm = zeros(L_segm, I); % Preallocate matrix to store DFT results for each segment

    % Compute DFT for each segment
    for i = 1:I
        segm_0 = (i-1) * L_segm + 1;
        segm_1 = segm_0 + L_segm - 1;
        x_segm = x(segm_0:segm_1); % Get ith segment
        dft_segm(:, i) = fft(x_segm); % Compute DFT of the ith segment (I no. of segments are distributed along columns)
    end
    
    dft_means = mean(dft_segm, 2); % Compute the average of each DFT value across all segments
    f_seg = (0:L_segm/2-1)*(Fs/L_segm); % Get frequency vector corresp. to each segment
    % Calculate FC for a set of discrete frequencies f1, f2, f3, f4
    % Numerator of FC
    mean_dev =      (dft_segm(f_seg == fv(1),:) - dft_means(f_seg == fv(1))) .* ...
                    (dft_segm(f_seg == fv(2),:) - dft_means(f_seg == fv(2))) .* ...
                    (dft_segm(f_seg == fv(3),:) - dft_means(f_seg == fv(3))) .* ...
                    (dft_segm(f_seg == fv(4),:) - dft_means(f_seg == fv(4))); % DFT deviations from the mean (in each segment)
    sum_mean_dev =  sum(mean_dev); % sum over all I segments
    fourth_pwr_mean_dev =  sum_mean_dev^4; % raise to the fourth power
    % Denominator of FC
    den_f1 = sum( (dft_segm(f_seg == fv(1),:) - dft_means(f_seg == fv(1))) .^4 );
    den_f2 = sum( (dft_segm(f_seg == fv(2),:) - dft_means(f_seg == fv(2))) .^4 );
    den_f3 = sum( (dft_segm(f_seg == fv(3),:) - dft_means(f_seg == fv(3))) .^4 );
    den_f4 = sum( (dft_segm(f_seg == fv(4),:) - dft_means(f_seg == fv(4))) .^4 );
    
    % Get FC
    FC = fourth_pwr_mean_dev / (den_f1*den_f2*den_f3*den_f4);
end
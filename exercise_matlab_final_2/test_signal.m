function [x, t] = test_signal(f0, Fs, L)

    % Ex:1
    % To create a digital test signal consisting of 4 random sine waves with frequencies f0,
    % 2f0, 3f0 and 4f0 (f0 is an arbitrary constant frequency); constant amplitudes of these
    % waves are: A1=1, A2=0.3, A3=0.2, A4=0.1 respectively and all waves have random initial
    % phases, uniformly distributed in range 0,2*pi
    
    % Define parameters
    A1 = 1.0;        % Amplitude of the first sine wave
    A2 = 0.3;        % Amplitude of the second sine wave
    A3 = 0.2;        % Amplitude of the third sine wave
    A4 = 0.1;        % Amplitude of the fourth sine wave
    
    % Generate random initial phases uniformly distributed between 0 and 2*pi
    phi1 = 2*pi*rand();
    phi2 = 2*pi*rand();
    phi3 = 2*pi*rand();
    phi4 = 2*pi*rand();
    
    % Generate the time vector
    T = 1/Fs;
    t = (0:L-1)*T;   % Time vector
    
    % Create the individual sine wave components
    sine_wave_1 = A1 * sin(2*pi*f0*t + phi1);
    sine_wave_2 = A2 * sin(2*pi*2*f0*t + phi2);
    sine_wave_3 = A3 * sin(2*pi*3*f0*t + phi3);
    sine_wave_4 = A4 * sin(2*pi*4*f0*t + phi4);
    
    % Sum up all the sine waves to create the final signal
    x = sine_wave_1 + sine_wave_2 + sine_wave_3 + sine_wave_4;

end
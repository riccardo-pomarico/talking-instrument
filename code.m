%% DAAP Homework #1 Group composition: 
% Marazzi Alice alice.marazzi@mail.polimi.it
% Pomarico Riccardo riccardo.pomarico@mail.polimi.it

close all
clearvars
clc

% Setting the number of taps

M_voice = 16;
M_music = 64;

% Loading audio files

[voice, ~] = audioread('speech.wav');
[music, Fs] = audioread('piano.wav');

% Defining the intervals

len = min(length(voice), length(music));
voice = voice(1:len);
music = music(1:len);

% Reducing peaks of the signals

music = 0.4 * music; 
voice = 0.4 * voice;

% Initializing output

y_closed = zeros(len,1);
y_steep = zeros(len,1);

% Defining the window variable

window_len = 1024;

% 50% of overlap

hop_size = window_len/2; 

blocks = floor((len - window_len) / hop_size) + 1;

% Hanning Window

win = hann(window_len);

for i = 1 : blocks-1

    % Windowing
    
    music_block = win.*music(i * hop_size+1 : i * hop_size+window_len);
    voice_block= win.*voice(i * hop_size+1 : i * hop_size+window_len);

    % FFT

    music_fft = fft(music_block);

    % Autocorrelation matrix and vector

    [autocorr_music, m_lags] = xcorr(music_block); 
    [autocorr_voice, v_lags] = xcorr(voice_block);

    % Cropping the vector

    p_music = autocorr_music(m_lags >= 1); 
    p_music = p_music(1:M_music);

    p_voice = autocorr_voice(v_lags >= 1); 
    p_voice = p_voice(1:M_voice);

    % Creating the autocorrelation vector needed for the Toeplitz matrix

    autocorr_music = autocorr_music(m_lags >= 0);
    autocorr_music = autocorr_music(1:M_music); 
    R_music = toeplitz(autocorr_music);

    autocorr_voice = autocorr_voice(v_lags >= 0); 
    autocorr_voice = autocorr_voice(1:M_voice);
    R_voice = toeplitz(autocorr_voice);

    % Closed form solutions

    w_music_closed = R_music \ p_music;
    w_voice_closed = R_voice \ p_voice;

    % Steepest Descent

    eig_m = eig(R_music);
    eig_v = eig(R_voice);

    % The mu value is computed to grant stability

    mu_m = 0.95*2/max(eig_m);
    mu_v = 0.5*2/max(eig_v);

    tau_m = 1/(2*mu_m*min(eig_m));
    tau_v = 1/(2*mu_v*min(eig_v));

    % Multiplying tau by 5 to reach the point of convergence

    iterations_m = 5 * tau_m;
    iterations_v = 5 * tau_v;

    % Limiting the number of iterations

    if (iterations_m>1000)
        iterations_m = 1000;
    end 

    if (iterations_v>1000)
        iterations_v = 1000;
    end 
    
    w_music_steep = zeros(size(p_music)); 
    w_voice_steep = zeros(size(p_voice));
     
    for n = 1:iterations_m
        grad = 2 * (p_music - R_music * w_music_steep(:, n));
        w_temp = w_music_steep(:, n) + 0.5 * mu_m * grad;
        w_music_steep = [w_music_steep, w_temp];
    end

    for n = 1:iterations_v
        grad = 2 * (p_voice - R_voice * w_voice_steep(:, n));
        w_temp = w_voice_steep(:, n) + 0.5 * mu_v * grad;
        w_voice_steep = [w_voice_steep, w_temp];
    end

    % Inverse filter coefficients

    A_m_closed = [1; -w_music_closed];
    A_m_closed_fft = fft(A_m_closed, size(music_block, 1));

    A_m_steep = [1; -w_music_steep(:, end)];
    A_m_steep_fft = fft(A_m_steep, size(music_block, 1));

    % Forward filter coefficients

    H_v_closed = [1; -w_voice_closed];
    H_v_closed_fft = fft (H_v_closed, size(voice_block, 1));
    
    H_v_steep = [1; -w_voice_steep(:, end)];
    H_v_steep_fft = fft (H_v_steep, size(voice_block, 1));

    % Whitening filter

    e_m_closed = music_fft.*A_m_closed_fft;
    e_m_steep = music_fft.*A_m_steep_fft;

    % Shaping filter

    out_closed_fft = e_m_closed./H_v_closed_fft;
    out_closed = ifft(out_closed_fft);
    out_steep_fft = e_m_steep./H_v_steep_fft;
    out_steep = ifft(out_steep_fft);

    % Computing overlap and add

    y_closed(i*hop_size + 1 : i*hop_size + window_len) = y_closed(i*hop_size + 1 : i*hop_size + window_len) + out_closed;
    y_steep(i*hop_size + 1 : i*hop_size + window_len) = y_steep(i*hop_size + 1 : i*hop_size + window_len) + out_steep;

end

% Cola condition

    % Initialize the sum of analysis windows

    sum_win = zeros(len, 1);
    
    % Loop over all frames and sum the analysis windows

    for i = 1:blocks
        idx = (i-1)*hop_size + (1:window_len);
        sum_win(idx) = sum_win(idx) + win;
    end
    
    % Check the COLA condition

    if round(sum_win(window_len-hop_size:end-window_len+hop_size))-(sum_win(window_len-hop_size:end-window_len+hop_size)) < 0.01
        disp('The window satisfies the COLA condition');
    else
        disp('The window does not satisfy the COLA condition');
    end

% Normalizing the output

y_closed = y_closed / (max(abs(y_closed)));
y_steep = y_steep / (max(abs(y_steep)));

audiowrite('Marazzi_Pomarico_result_closedform.wav', y_closed, Fs);
audiowrite('Marazzi_Pomarico_result_steepestdescent.wav', y_steep, Fs);
function rho = stgi(clean_speech, degraded_speech, Fs)
%     rho = stgi(clean_speech, degraded_speech, Fs) returns the output of 
%     the Spectro-Temporal Glimpsing Index (STGI) predictor.
% 
%     Implementation of the Spectro-Temporal Glimpsing Index (STGI) 
%     predictor, described in [1].
% 
%     Inputs:
%            clean_speech:        clean reference time domain signal
%            degraded_speech:     noisy/processed time domain signal
%            Fs:                  sampling rate [Hz]
% 
%     Output:
%            rho: intelligibility index

%     References:
%     [1] A. Edraki, W.-Y. Chan, J. Jensen, & D. Fogerty, 
%     “A Spectro-Temporal Glimpsing Index (STGI) for Speech Intelligibility
%     Prediction," Proc. Interspeech, 5 pages, Aug 2021.


%     Copyright 2021, Amin Edraki, a.edraki@queensu.ca, Queen's University
%     at Kingston, Multimedia Coding and Communications Laboratory.

%     stgi is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     stgi is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with wstmi.  If not, see <https://www.gnu.org/licenses/>.

if length(clean_speech) ~= length(degraded_speech)
    error('clean and degraded signals should have the same length!');
end

% resampling the input signals
fs = 10000;
if Fs ~= fs
    clean_speech = resample(clean_speech, fs, Fs);
    degraded_speech = resample(degraded_speech, fs, Fs);
end


STM_channels = ones(11, 4);
thresholds   = [0.252,0.347,0.275,0.189;
    0.502,0.495,0.404,0.279;
    0.486,0.444,0.357,0.247;
    0.456,0.405,0.332,0.229;
    0.426,0.361,0.287,0.191;
    0.357,0.299,0.229,0.150;
    0.269,0.228,0.175,0.114;
    0.185,0.158,0.118,0.075;
    0.119,0.103,0.073,0.047;
    0.081,0.067,0.047,0.030;
    0.050,0.043,0.031,0.020];

% Parameters of the log Mel-spectrogram
win_length  = 25.6;
win_shift   = win_length/2;
freq_range  = [64 min(floor(fs./2), 12000)];
num_bands   = 130;
band_factor = [];
N           = 40;

% Parameters of the voice activity detector
dyn_range   = 40;
N_frame     = 256;



% Voice activity detection
[clean_speech, degraded_speech] = removeSilentFrames(clean_speech, degraded_speech, dyn_range, N_frame, N_frame/2);


% log Mel-spectrograms
[X, ~]    = log_mel_spectrogram(clean_speech, fs, win_shift, win_length, freq_range, num_bands, band_factor);
[Y, ~]    = log_mel_spectrogram(degraded_speech, fs, win_shift, win_length, freq_range, num_bands, band_factor);

% Spectro-temporal modulation decomposition
X_hat     = sgbfb(X, STM_channels);
Y_hat     = sgbfb(Y, STM_channels);




% Initializing the intelligibility vector
intell    = zeros(length(N:size(X_hat, 4)), 1);

for n = N:size(X_hat, 4)

    X_seg = X_hat(:, :, :, (n-N+1):n);
    Y_seg = Y_hat(:, :, :, (n-N+1):n);
    
    % Row normalization
    X_seg = X_seg - mean(X_seg, 4);
    Y_seg = Y_seg - mean(Y_seg, 4);
    X_seg = X_seg ./ sqrt(sum(X_seg.*X_seg, 4));
    Y_seg = Y_seg ./ sqrt(sum(Y_seg.*Y_seg, 4));
    
    % Column normalization
    X_seg = X_seg - mean(X_seg, 3);
    Y_seg = Y_seg - mean(Y_seg, 3);
    X_seg = X_seg ./ sqrt(sum(X_seg.*X_seg, 3));
    Y_seg = Y_seg ./ sqrt(sum(Y_seg.*Y_seg, 3));
     
    % Similarity measure
    d = squeeze(sum(X_seg .* Y_seg, 3));
    d = squeeze(mean(d, 3));
    
    % Glimpse detection
    g = d > thresholds;
    intell(n-N+1) = mean(g(:));
    

end
rho = mean(intell(:));
end







function [x_sil y_sil] = removeSilentFrames(x, y, range, N, K)
%     [X_SIL Y_SIL] = REMOVESILENTFRAMES(X, Y, RANGE, N, K) X and Y
%     are segmented with frame-length N and overlap K, where the maximum energy
%     of all frames of X is determined, say X_MAX. X_SIL and Y_SIL are the
%     reconstructed signals, excluding the frames, where the energy of a frame
%     of X is smaller than X_MAX-RANGE

%     Copyright 2009: Delft University of Technology, Signal & Information
%     Processing Lab. The software is free for non-commercial use. This program
%     comes WITHOUT ANY WARRANTY.

x       = x(:);
y       = y(:);

frames  = 1:K:(length(x)-N);
w       = hanning(N);
msk     = zeros(size(frames));

for j = 1:length(frames)
    jj      = frames(j):(frames(j)+N-1);
    msk(j) 	= 20*log10(norm(x(jj).*w)./sqrt(N));
end

msk     = (msk-max(msk)+range)>0;
count   = 1;

x_sil   = zeros(size(x));
y_sil   = zeros(size(y));

for j = 1:length(frames)
    if msk(j)
        jj_i            = frames(j):(frames(j)+N-1);
        jj_o            = frames(count):(frames(count)+N-1);
        x_sil(jj_o)     = x_sil(jj_o) + x(jj_i).*w;
        y_sil(jj_o)  	= y_sil(jj_o) + y(jj_i).*w;
        count           = count+1;
    end
end

x_sil = x_sil(1:jj_o(end));
y_sil = y_sil(1:jj_o(end));
end




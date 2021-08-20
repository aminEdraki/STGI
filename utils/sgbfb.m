function STM_dec = sgbfb(log_mel_spec, STM_channels)
%     Separated spectro-temporal modulation Gabor filter bank for stgi.
%     This script extracts the separable spectro-temporal modulation
%     represntation from a log Mel-spectrogram as described in [1, 2].
%     Modified from the reference [1, 3] implementation.

%     Inputs:
%       log_mel_spec    log Mel-spectrogram
%       size_max        maximum filter size (lower center modulation 
%       frequency) in [bands, frames]
%       nu              half-waves under the envelope in [spectral 
%       temporal] dimension
%       phases          Phase of [spectral temporal] modulation filters
%       STM_channels    A binary matrix selecting the STM channels to be
%       calculated

%     Output:
%       STM_dec:        spectro-temporal modulation decomposition output


%     References:
%     [1] Schädler, M. R., & Kollmeier, B. (2015). Separable spectro-
%     temporal Gabor filter bank features: Reducing the complexity of 
%     robust features for automatic speech recognition. The Journal of the  
%     Acoustical Society of America, 137(4), 2047-2059.
%     [2] Edraki, A., Chan, W. Y., Jensen, J., & Fogerty, D. (2020). 
%     Speech Intelligibility Prediction Using Spectro-Temporal Modulation
%     Analysis. IEEE/ACM Transactions on Audio, Speech, and Language 
%     Processing, 29, 210-225.
%     [3] Schädler, M. R., Meyer, B. T., & Kollmeier, B. (2012). 
%     Spectro-temporal modulation subspace-spanning filter bank features 
%     for robust automatic speech recognition. The Journal of the 
%     Acoustical Society of America, 131(5), 4134-4151.




%     Copyright (C) 2015-2018 Marc RenÃ© SchÃ¤dler
%     E-mail marc.r.schaedler@uni-oldenburg.de
%     Institute Carl-von-Ossietzky University Oldenburg, Germany
% 
%     Copyright 2021, Amin Edraki, a.edraki@queensu.ca, Queen's University
%     at Kingston, Multimedia Coding and Communications Laboratory.

%     This file is part of wstmi.
% 
%     wstmi is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     wstmi is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with wstmi.  If not, see <https://www.gnu.org/licenses/>.

%% Default settings and checks

% Max spectral and temporal modulation frequencies
omega_max = [2*pi/3 pi/2];

% Get number of bands from input
num_bands = size(log_mel_spec,1);

% Default number of half-waves under the envelope
nu = [3.5 3.5];


% Default maximum filter size (lower center modulation frequencies)
size_max = [3*num_bands 40];

% Default Phases of modulation filters [spectral temporal]
phases = [0 0];

% Distance between the adjacent filters
distance = [0.2 0.2];

% Default spectro temporal channels
if isempty(STM_channels)
  STM_channels = ones(11, 5);
end

% Maximum context
context = floor(size_max(2)/2);




  


%% Calculate separated Gabor filter bank features

  
% Note: The order of spectral and temporal modulation processing
% does not matter

[row, col] = find(STM_channels);
row = unique(row);
col = unique(col);

% Calculating modulation center frequencies
omega_s = gbfb_axis(omega_max(1),size_max(1),nu(1),distance(1));
omega_t = gbfb_axis(omega_max(2),size_max(2),nu(2),distance(2));

% Including the DC filter
omega_s = [0, omega_s];
omega_t = [0, omega_t];

% Selecting the channels
omega_s = omega_s(row);
omega_t = omega_t(col);


% Temporally pad log Mel-spectrogram by repeating first and last frames
log_mel_spec  = [repmat(log_mel_spec(:,1),1,context) log_mel_spec repmat(log_mel_spec(:,end),1,context)];

% Spectral 1D filtering
features_spec = gbfb1d(log_mel_spec, size_max(1), nu(1), phases(1), omega_s);
% Temporal 1D filtering
features_spec = gbfb1d(features_spec.', size_max(2), nu(2), phases(2), omega_t).';


% Initializing the spectro-temporal decomposition matrix
time_frames   = size(log_mel_spec, 2);
freq_bins     = size(log_mel_spec, 1);
STM_dec       = zeros(length(omega_s), length(omega_t), freq_bins, time_frames-2*context);

for s = 1:length(omega_s)
    for r = 1:length(omega_t)
        filtered_spec = features_spec((s-1)*freq_bins+1:s*freq_bins, (r-1)*time_frames+1:r*time_frames);
        STM_dec(s, r, :, :) = filtered_spec(:, (1+context):(end-context));
    end
end

end
 

%%%%%%%%%% 1D Gabor filter bank (1D-GBFB) functions %%%%%%%%%%
% Most of this code is modified from the reference GBFB implementation [2]

function out = gbfb1d(in, size_max, nu, phase, omega)
% One dimensional Gabor filter bank filtering

% For each center modulation frequency
out = [];
for i=1:length(omega)
  % Generate a 1D Gabor filter
  gfilter = gfilter_gen(omega(i), nu, phase, size_max);
  % Filter the input
  in_filtered = conv2(in, gfilter, 'same');
  
  out = [out; in_filtered];
end
end



function gfilter = gfilter_gen(omega, nu, phi, size_max)
% Generates a 1D Gabor filter function
% omega       Center frequency in rad
% nu          Number of half-waves under the envelope
% phi         Phase of gfilter at the center sample in rad
% size_max    Maximum size to determine when to use the envelope
%             as the filter function

w = inf;
if omega > 0
    w = 2*pi / abs(omega) * nu / 2;
end
if w > size_max
    w = size_max;
    omega = 0;
end
envelope = hann_win(w); % C.f. Equation 1a in [1]
win_size = length(envelope);
x_0 = (win_size+1) / 2;
x = 1:win_size;
sinusoid = exp(1i.*(omega*(x - x_0) + phi)); % C.f. Equation 1b in [1]
gfilter  = real(envelope(:) .* sinusoid(:)); % C.f. Equation 1c in [1]
envelope_mean = mean(mean(envelope));
gfilter_mean = mean(mean(gfilter));
if (omega ~= 0)
    gfilter = gfilter - envelope./envelope_mean .* gfilter_mean;
end
gfilter = gfilter ./ max(abs(fft(gfilter)));


end



function window_function = hann_win(width)
% A hanning window of "width" with the maximum centered on the center sample
x_center = 0.5;
step = 1/width;
right = x_center:step:1;
left = x_center:-step:0;
x_values = [left(end:-1:1) right(2:end)].';
valid_values_mask = (x_values > 0) & (x_values < 1);
window_function = 0.5 * (1 - ( cos(2*pi*x_values(valid_values_mask))));
end


function omega = gbfb_axis(omega_max, size_max, nu, distance)
% Calculates the center modulation frequencies (c.f. [2])
% omega_max   Maximum center modulation frequency in rad
% size_max    Maximum extension of a Gabor filter in samples
% nu          Number of half-waves under the envelope
% distance    Distance between adjacent filters
omega_min = (pi * nu) / size_max;
c = distance * 8 / nu;
space = (1 + c/2) / (1 - c/2);
count = 0;
omega(1) = omega_max;
while omega(end)/space > omega_min
  omega(1+count) = omega_max/space^count;
  count = count + 1;
end
omega = fliplr(omega);
end




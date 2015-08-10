% demonstration of the waveform train decomposition
clc; clear all; close all;

disp 'reading EEG ...'

% Load Kari's clips
%{
filename = 'E:\data\human CNS\EMD\Sz\clips\2012PP05Sz1.mat';
signal = load(filename);
signal = signal.data;
[numChannels, numSamples] = size(signal); % each segment is 10 minutes
fs = numSamples/(10*60);
channelsToUse = [1:4 6:24 26:31 35:52 54:59 61:87 89:94 96:98];   % selected a subset of channels
signal = signal(channelsToUse,:);
%}

% Load from EDF
%
filename = 'data/KT_7.edf';
hdr = edfopen(filename);
signal = edfread(hdr,0,hdr.nsamples);
channelsToUse = [2:5 7:25];   % selected a subset of channels
fs = hdr.samples_per_second;
signal = signal(channelsToUse,:);
%}

disp 'filtering...'
% notch filter
[b,a] = iirnotch(60/(fs/2),0.5/fs); % determie filter coefficients
signal = filtfilt(b,a,signal')'; % apply filter

% high-pass filter
cutoff = 2; % Hz
k = hamming(round(fs/cutoff)*2+1); % hamming window
k = k/sum(k); % normalize
signal = signal - convmirr(signal',k)'; % filter

% Doesn't work without hdr (not available from Kari's clips). 
%
disp 'plotting raw data...'
t = (0:size(signal,2)-1)/fs;
yticks = 1e4*(1:size(signal,1));
plot(t,bsxfun(@plus,signal',yticks))
set(gca,'YTick',yticks,'YTickLabel',arrayfun(@(i) strtrim(hdr.channelnames(i,:)),channelsToUse, 'uni', false))
xlabel 'time (s)'
%}

% algorithm paramaters (all units are in samples)
startTime =270; % (s)
epoch = 2500;  % samples
epochStep = 2000;
waveform_width = 101;
ntrains = 3;

iterations = round(startTime*fs):epochStep:size(signal,2)-epoch;
waveforms = zeros(waveform_width, length(channelsToUse),ntrains,length(iterations)); % stores all waveforms
occurrences = zeros(epoch, ntrains, length(iterations)); % stores all occurrences

k = 1;
for i=iterations
    segment = signal(:,i+(1:epoch))'; % segment is the raw data
    [w, u] = choo3(segment, ntrains, waveform_width); % u = occurrences, w = waveforms.
    waveforms(:,:,:,k) = w;
    occurrences(:,:,k) = u;
    %show_trains(segment, u, w)
    %keyboard
    disp([num2str(k), ' of ', num2str(length(iterations))]);
    k = k + 1;
end

%% Example of how to plot saved data
show_iteration = 3;
show_trains(signal(:,iterations(show_iteration)+(1:epoch))', ... 
    occurrences(:,:,show_iteration), waveforms(:,:,:,show_iteration))

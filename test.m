% demonstration of the waveform train decomposition

disp 'reading EEG ...'
filename = 'data/KT_7.edf';
hdr = edfopen(filename);
signal = edfread(hdr,0,hdr.nsamples);
channelsToUse = [5];
% channelsToUse = [2:5 7:25];   % selected a subset of channels
% channelsToUse = [2:5 7:9];
fs = hdr.samples_per_second;
signal = signal(channelsToUse,:);

for i = 1:2
    signal(i,:) = signal(1,:);
end

disp 'filtering...'
% notch filter
[b,a] = iirnotch(60/(fs/2),0.5/fs); % determie filter coefficients
signal = filtfilt(b,a,signal')'; % apply filter

% high-pass filter
cutoff = 2; % Hz
k = hamming(round(fs/cutoff)*2+1); % hamming window
k = k/sum(k); % normalize
signal = signal - convmirr(signal',k)'; % filter

disp 'plotting raw data...'
t = (0:size(signal,2)-1)/fs;
yticks = 2e4*(1:size(signal,1));
plot(t,bsxfun(@plus,signal',yticks))
set(gca,'YTick',yticks,'YTickLabel',arrayfun(@(i) strtrim(hdr.channelnames(i,:)),channelsToUse, 'uni', false))
xlabel 'time (s)'

% algorithm paramaters (all units are in samples)
startTime = 270; % (s)
epoch = 2500;  % samples
epochStep = 2000;
waveform_width = 101; % must be odd
ntrains = 1;

numTrains = 5;

for i=round(startTime*fs):epochStep:size(signal,2)-epoch
    segment = signal(:,i+(1:epoch))'; % segment is the raw data
    segment_copy = segment;
    [T, number_of_channels] = size(segment);  % T = # time samples
    u_all = zeros(T, numTrains);    % occurrences
    w_all = zeros(waveform_width,number_of_channels, numTrains);  % waveforms
    for j= 1:numTrains
        [w, u] = choo3(segment_copy, ntrains, waveform_width);
        u_all(:,j) = squeeze(u);
        w_all(:,:,j) = squeeze(w);
        show_trains(segment_copy(:,1), u, w(:,1,:))
        y = reconstruct(u,w);
        segment_copy = segment_copy - y;
    end
    show_trains(segment(:,1), u_all, w_all(:,1,:));
end
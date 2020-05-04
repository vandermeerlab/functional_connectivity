function [X, cfg_out] = amp_crosscorr_on_gamma_events_new(cfg_in, CSC1, CSC2, eventStats, eventData, eventName, varargin)
% original code from Adhikari,Sigurdsson,Topiwala,& Gordon, 2010. J.NsciMethods
% 2014-09-24. JJS. Modified so that it can run through multiple sessions
% and save the data and figures to each session folder.
% 2020-04-17. JJS. Modified to work within vdmlab codebase.
% Operates on a single session.

% INPUTS:
%   cfg_in:                         config input
%       cfg_in.cfg_out.low_freq             low cut off, in Hz, of the band pass filter that will be applied to CSC1.data and CSC2.data
%       cfg_in.cfg_out.high_freq            high cut off, in Hz, of the band pass filter that will be applied to CSC1.data and CSC2.data
%   CSC1                            tsd containing local field potential from brain area 1
%   CSC2                            tsd containing local field potential from brain area 2
%   eventStats:                     IV structure with start and end timestamps of LFP events. Also contains information about eventData.
%       evenStats.fs                    sampling rate of the data

% OUTPUTS:
%   lags:                           vector contaning lags from -100 ms to +100 ms, over which the crosscorrelation was done
%   crosscorr:                      vector with the crosscorrelation of the amplitude of CSC1.data CSC2.data after being filtered between cfg_out.low_freq and cfg_out.high_freq
%   max_crosscorr_lag:              lag at which the crosscorrelation peaks. Negative max_crosscorr_lag indicates that CSC1.data is leading CSC2.data.
%   g:                              index of the max crosscorrelation lag
%   cfg_out:

doPlotKeeps = 1;
doPlotEverything = 0;
Inspect = 1;
doSave = 1;
process_varargin(varargin);

cfg_def = [];
samp_freq = eventStats.usr.fs;
eventDuration = eventStats.cfg.history.cfg.eventDuration;
dt = 1/samp_freq;
cfg_def.percentile = 0.95;
cfg_def.nShuffle = 1000;
cfg_def.maxlags = round(eventDuration/dt);              % number of lags to calculate on either side of zero. Duration of maxlags = 1/2 LFP eventDuration.
cfg_def.low_freq = 45;
cfg_def.high_freq = 65;
cfg_out = ProcessConfig(cfg_def, cfg_in);

%% Checks 
if length(CSC1.data)~= length(CSC2.data);
    error('ERROR in amp_crosscorr. CSC1.data and CSC2.data must be vectors of the same size;')
end
s = size(CSC1.data);
if min(s) ~= 1
    error('ERROR in amp_crosscorr. CSC1.data and CSC2.data must be one-dimensional vectors')
end
s = size(CSC2.data);
if min(s) ~= 1
    error('ERROR in amp_crosscorr. CSC1.data and CSC2.data must be one-dimensional vectors')
end

%% Calculate the amplitdue envelope
order = round(samp_freq); % determines the order of the filter used
if mod(order, 2) ~= 0
    order = order-1;
end
Nyquist = floor(samp_freq/2); % determines nyquist frequency
MyFilt = fir1(order, [cfg_out.low_freq cfg_out.high_freq]/Nyquist); %creates filter
filtered1 = Filter0(MyFilt, CSC1.data); % filters CSC1.data between cfg_out.low_freq and cfg_out.high_freq
filtered2 = Filter0(MyFilt, CSC2.data); % filters CSC2.data between cfg_out.low_freq and cfg_out.high_freq
filt_hilb1 = hilbert(filtered1); % calculates the Hilbert transform of CSC1.data
amp1 = abs(filt_hilb1); % calculates the instantaneous amplitude of CSC1.data filtered between cfg_out.low_freq and cfg_out.high_freq
amp1 = amp1 - mean(amp1); % removes mean of the signal because the DC component of a signal does not change the correlation
filt_hilb2 = hilbert(filtered2); % calculates the Hilbert transform of CSC2.data
amp2 = abs(filt_hilb2); % calculates the instantaneous amplitude of CSC2.data filtered between cfg_out.low_freq and cfg_out.high_freq
amp2 = amp2 - mean(amp2);
assert(length(eventStats.usr.istart) == length(eventStats.usr.iend))

iEvent = 0;
tic   % takes about 5 minutes for 1000 LFP events (April 2020)
for iEvent = 1:length(eventStats.usr.istart);
    disp(strcat(num2str(iEvent), ' of_', num2str(length(eventStats.usr.istart))))    
    amp1touse = amp1(eventStats.usr.istart(iEvent): eventStats.usr.iend(iEvent));
    amp2touse = amp2(eventStats.usr.istart(iEvent): eventStats.usr.iend(iEvent));
    [crosscorr(:,iEvent), lags] = xcorr(amp1touse, amp2touse, cfg_out.maxlags, 'coeff'); % calculates crosscorrelations between amplitude vectors
    lags = (lags./samp_freq)*1000; % converts lags to miliseconds
    [C(iEvent), g(iEvent)] = max(crosscorr(:,iEvent));
    max_crosscorr_lag(iEvent) = lags(g(iEvent)); % identifies the lag at which the crosscorrelation peaks
    
    parfor iShuf = 1: cfg_out.nShuffle
        amp2_shuff = AUX_shuffle_phases_new(amp2touse');
        temp_corrvalues = xcorr(amp1touse, amp2_shuff, cfg_out.maxlags, 'coeff');
        shuff_max_xcorr(iShuf) = max(temp_corrvalues);
    end
    sorted_shuff = sort(shuff_max_xcorr);
    percentile_index = cfg_out.percentile*length(sorted_shuff);
    max_xcorr_cutoff(iEvent) = sorted_shuff(round(percentile_index));
    keep(iEvent) = C(iEvent) > max_xcorr_cutoff(iEvent);
end
toc

X.lags = lags;
X.crosscorr = crosscorr;
X.max_crosscorr_lag = max_crosscorr_lag;
X.C = C;
X.g = g;
X.max_xcorr_cutoff = max_xcorr_cutoff;
X.keep = keep;
X.cfg = cfg_out;
X.eventData = eventData;    % add the actual data for the LFP events to the output X, in case needed later
X.eventStats = eventStats;  % add the intervals for LFP events to the output X, in case needed later

if doSave == 1;
    SSN = GetSSN('SingleSession');
    fn = strcat(SSN, '-XCorr-', eventName);
    save(fn, 'X');
    disp(fn);
    disp('data saved');
end

if doPlotKeeps == 1;
    figure(1); clf
    for iEvent = 1:length(eventStats.usr.istart);
        plot(X.lags, X.crosscorr(:,iEvent),'color',[0 0 1],'linewidth',2),hold on % plots crosscorrelations
        plot(X.lags(X.g(iEvent)), X.crosscorr(X.g(iEvent), iEvent),'rp','markerfacecolor',[1 0 0],'markersize',10) % plots marker at the peak of the cross correlation
        plot([0 0],[1.05*max(X.crosscorr(:,iEvent)) 0.95*min(X.crosscorr(:,iEvent))],'color',[0 0 0],'linestyle',':', 'linewidth',2) % plots dashed line at zero lag
        set(gca,'xtick',[-100 -50 0 50 100])
        axis tight, box off, xlim([-101 100])
        xlabel('Lag (ms)','fontsize',14)
        ylabel('Crosscorrelation','fontsize',14)
        %             disp('press any key to continue')
        title(num2str(iEvent))
        if Inspect
            pause
            clf
        end
    end
end

if doPlotEverything == 0;
%     clf
%     % crosscorr plot fof this gamma event
%     subplot(2,2,1)
%     plot(lags, crosscorr,'color',[0 0 1],'linewidth',2),hold on %plots crosscorrelations
%     plot(lags(g),crosscorr(g),'rp','markerfacecolor',[1 0 0],'markersize',10)%plots marker at the peak of the cross correlation
%     plot([0 0],[1.05*max(crosscorr) 0.95*min(crosscorr)],'color',[0 0 0],'linestyle',':', 'linewidth',2) %plots dashed line at zero lag
%     set(gca,'xtick',[-100 -50 0 50 100])
%     axis tight, box off, xlim([-101 100])
%     xlabel('Lag (ms)','fontsize',14)
%     ylabel('Crosscorrelation','fontsize',14)
%     title(strcat('max crosscorr = ', num2str(max_crosscorr_value_real)))
%     % crosscorr values for vstr time series vs. 1,000 phase permuted ofc time series
%     subplot(2,2,2)
%     hist(max_crosscorr_value_shuff, 50);
%     c = axis;
%     line([max_crosscorr_value_real max_crosscorr_value_real], [c(3) c(4)], 'Color', 'r', 'LineWidth', 1);
%     title(strcat('keep = ', num2str(keep(iEvent))))
%     xlabel('max crosscorr value', 'fontsize', 14)
%     ylabel('count', 'fontsize', 14)
%     % raw LFP traces for OFC and VSTR
%     subplot(2,2,3)
%     CSC1.data = X(1,:,iEvent);
%     CSC2.data = X(2,:,iEvent);
%     hold on
%     plot(CSC1.data, 'r')  % ofc
%     plot(CSC2.data, 'b')  % vstr
%     legend('ofc', 'vstr')
%     % g50 amplitude envelopes for OFC and VSTR
%     subplot(2,2,4)
%     Gofc = tsd(lags(1:end-1)', CSC1.data');
%     Gvstr = tsd(lags(1:end-1)', CSC2.data');
%     [~, ~, ~, ~, IEofc]=InstSig(Gofc,45,65,64);
%     [~, ~, ~, ~, IEvstr]=InstSig(Gvstr,45,65,256);
%     plot(Gofc.range, IEofc.data, 'r')  % ofc
%     plot(Gvstr.range, IEvstr.data, 'b')  % vstr
%     
%     pause
end

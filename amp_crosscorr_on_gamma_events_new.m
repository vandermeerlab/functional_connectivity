function [lags, crosscorr, max_crosscorr_lag, g, C, cfg_out]=amp_crosscorr_on_gamma_events_new(cfg_in, CSC1, CSC2, eventStats, varargin)
% original code from Adhikari,Sigurdsson,Topiwala,& Gordon, 2010. J.NsciMethods
% 2014-09-24. JJS. Modified so that it can run through multiple sessions
% and save the data and figures to each session folder.
% 2020-04-17. JJS. Modified to work within vdmlab codebase.
% Operates on a single session.

% INPUTS:
%   cfg_in:                         config input
%       cfg_in.cfg_out.low_freq             low cut off, in Hz, of the band pass filter that will be applied to eeg1 and eeg2
%       cfg_in.cfg_out.high_freq            high cut off, in Hz, of the band pass filter that will be applied to eeg1 and eeg2
%   CSC1                            tsd containing local field potential from brain area 1
%   CSC2                            tsd containing local field potential from brain area 2
%   eventStats:                     IV structure with start and end timestamps of LFP events. Also contains information about eventData.
%       evenStats.fs                    sampling rate of the data

% OUTPUTS:
%   lags:                           vector contaning lags from -100 ms to +100 ms, over which the crosscorrelation was done
%   crosscorr:                      vector with the crosscorrelation of the amplitude of eeg1 eeg2 after being filtered between cfg_out.low_freq and cfg_out.high_freq
%   max_crosscorr_lag:              lag at which the crosscorrelation peaks. Negative max_crosscorr_lag indicates that eeg1 is leading eeg2.
%   g:                              index of the max crosscorrelation lag
%   cfg_out:

doPlot = 1;
Inspect = 1;
doSave = 1; 
process_varargin(varargin);

cfg_def = [];
samp_freq = eventStats.usr.fs;
eventDuration = eventStats.cfg.history.cfg.eventDuration;
dt = 1/samp_freq;
cfg_def.percentile = 0.95;
cfg_def.nShuffle = 1000;
cfg_def.maxlags = round(eventDuration/dt);              % number of lags to calculate on either side of zero. Duration of maxlags = LFP eventDuration.
cfg_def.low_freq = 45;
cfg_def.high_freq = 65;
cfg_out = ProcessConfig(cfg_def, cfg_in);

fd = FindFiles('*keys.m'); max_crosscorr_lag = nan(1,length(fd));
for iSess = 1:length(fd);
    pushdir(fileparts(fd{iSess}));
    SSN = GetSSN('SingleSession'); disp(SSN);
    
    eeg1 = CSC1.data; eeg2 = CSC2.data; % samp_freq = 1/median(diff(CSC1.tvec));
    if length(eeg1)~= length(eeg2);
        error('ERROR in amp_crosscorr. eeg1 and eeg2 must be vectors of the same size;')
    end
    s=size(eeg1);
    if min(s)~=1
        error('ERROR in amp_crosscorr. eeg1 and eeg2 must be one-dimensional vectors')
    end
    s=size(eeg2);
    if min(s)~=1
        error('ERROR in amp_crosscorr. eeg1 and eeg2 must be one-dimensional vectors')
    end
    order = round(samp_freq); % determines the order of the filter used
    if mod(order,2)~= 0
        order = order-1;
    end
    Nyquist=floor(samp_freq/2); % determines nyquist frequency
    MyFilt=fir1(order,[cfg_out.low_freq cfg_out.high_freq]/Nyquist); %creates filter
    filtered1 = Filter0(MyFilt,eeg1); % filters eeg1 between cfg_out.low_freq and cfg_out.high_freq
    filtered2 = Filter0(MyFilt,eeg2); % filters eeg2 between cfg_out.low_freq and cfg_out.high_freq
    filt_hilb1 = hilbert(filtered1); % calculates the Hilbert transform of eeg1
    amp1 = abs(filt_hilb1); % calculates the instantaneous amplitude of eeg1 filtered between cfg_out.low_freq and cfg_out.high_freq
    amp1=amp1-mean(amp1); % removes mean of the signal because the DC component of a signal does not change the correlation
    filt_hilb2 = hilbert(filtered2); % calculates the Hilbert transform of eeg2
    amp2 = abs(filt_hilb2); % calculates the instantaneous amplitude of eeg2 filtered between cfg_out.low_freq and cfg_out.high_freq
    amp2=amp2-mean(amp2);
    assert(length(eventStats.usr.istart) == length(eventStats.usr.iend))
    
    tic
    for iEvent = 1:length(eventStats.usr.istart);
        disp(num2str(iEvent))
        amp1touse = amp1(eventStats.usr.istart(iEvent): eventStats.usr.iend(iEvent));
        amp2touse = amp2(eventStats.usr.istart(iEvent): eventStats.usr.iend(iEvent));
        [crosscorr(:,iEvent), lags] = xcorr(amp1touse, amp2touse, cfg_out.maxlags, 'coeff'); % calculates crosscorrelations between amplitude vectors
        lags = (lags./samp_freq)*1000; % converts lags to miliseconds
%         g(iEvent) = find(crosscorr(:,iEvent) == max(crosscorr(:,iEvent))); % identifies index where the crosscorrelation peaks
        [C(iEvent), g(iEvent)] = max(crosscorr(:,iEvent)); 
        max_crosscorr_lag(iEvent) = lags(g(iEvent)); % identifies the lag at which the crosscorrelation peaks
        
        if cfg_out.nShuffle > 0
            for iShuf = cfg_out.nShuffle:-1:1
                amp2_shuff = AUX_shuffle_phases_new(amp2touse');
                temp_corrvalues = xcorr(amp1touse, amp2_shuff, cfg_out.maxlags, 'coeff');
                shuff_max_xcorr(iShuf) = max(temp_corrvalues);
            end
        end
        sorted_shuff = sort(shuff_max_xcorr);
        percentile_index = cfg_out.percentile*length(sorted_shuff);
        max_xcorr_cutoff(iEvent) = sorted_shuff(round(percentile_index));
        eventsToUse(iEvent) = C(iEvent) > max_xcorr_cutoff(iEvent); 
        
        if doPlot == 1;
            figure(1); clf
            plot(lags, crosscorr(:,iEvent),'color',[0 0 1],'linewidth',2),hold on % plots crosscorrelations
            plot(lags(g(iEvent)),crosscorr(g(iEvent),iEvent),'rp','markerfacecolor',[1 0 0],'markersize',10) % plots marker at the peak of the cross correlation
            plot([0 0],[1.05*max(crosscorr(:,iEvent)) 0.95*min(crosscorr(:,iEvent))],'color',[0 0 0],'linestyle',':', 'linewidth',2) % plots dashed line at zero lag
            set(gca,'xtick',[-100 -50 0 50 100])
            axis tight, box off, xlim([-101 100])
            xlabel('Lag (ms)','fontsize',14)
            ylabel('Crosscorrelation','fontsize',14)
            %             disp('press any key to continue')
            title(num2str(iEvent))
            if Inspect
                pause
            end
        end    
    end
    toc
end

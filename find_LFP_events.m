function [eventStats, eventData] = find_LFP_events(cfg_in, fd, eventName)
% 2020-04-14. JJS. This function gets the start and end times of threshold crossing events for the filtered gamma LFP.

% Inputs:
%   cfg_in:                         config input
%   fd:                             list of sessions to analyze
%   eventName                       name (string) that the file will be saved as; for instance 'gamma50'

% Outputs:
%   eventStats:                     IV structure with event start and end times

%   eventData:                      nChannels x Samples x nEvents
% Typically, channel 1 = OFC. channel 2 = vStr     % can add Hipp later as Channel 3
% nSamples = number of data points in an LFP event (proportional to eventDuration)
% nEvents = number of events in the session

% Important Parameters
%   cfg_def.lowpass:                lowpass filter value (45 for g50, 70 for gamma50)
%   cfg_def.highpass:               highpass filter value (65 for g50, 90 for gamma80)
%   cfg_def.eventDuration:          how long you want the LFP events to be (in s). This parameter is chosen to make them all the same length, but in reality, evnets (for insteance, gamma events) will vary in duration.
%   cfg_def.percentile              the percentile value of LFP voltages above which you want to use as threshold for event detection


% fd = FindFiles('*keys.m');
% startSess = 1; endSess = length(fd);

cfg_def = [];
cfg_def.doPlot = 1;
cfg_def.MarkerSize = 20;
cfg_def.SaveIt = 0;
cfg_def.lowpass = 45;
cfg_def.highpass = 65;
cfg_def.eventDuration = 0.1;    % Duration (in seconds) of the gamma events. This is somewhat arbitrary. Julien used 0.4 s, Eric used 0.1 s
cfg_def.percentile = 0.95;
cfg_def.minIEI = .1;            % Duration (in seconds) that events need to be separated by.  IEI = inter-event-interval
cfg_out = ProcessConfig(cfg_def, cfg_in);

% for iSess = startSess:endSess;
for iSess = 1:length(fd);
    
    pushdir(fd{iSess});
    SSN = GetSSN('SingleSession'); disp(SSN);
    
    cfg_prep = [];
    cfg_prep.decimateByFactor = 2; % Downsampling the data. Default here is 2 (from 2kHz to 1kHz).
    cfg_prep.detrend = 1;   % for removing slow DC shifts in voltage
    cfg_prep.diffdata = 0;  % for removing autocorrelation in the time series
    cfg_prep.doRestrict = 1; % Restirct CSC to track time only.
    cfg_prep.hippflag = 0;  % Do we want to Load a hippocampal CSC? (for sessions with electrode in fissure). 0 = no. 1 = yes.
    [CSC_ofc, CSC_vstr, ~, ~] = prepCSCs_new(cfg_prep);   % can change this later to return hippocampus LFP also
    fs = 1/median(diff(CSC_vstr.tvec));  % sampling rate of the LFP
    
    %% Filter in the frequency band of interest
    cfg_filter = [];
    cfg_filter.f = [cfg_out.lowpass cfg_out.highpass];   % frequency band to filter in. So far, just using low gamma events [45 65]
    cfg_filter.display_filter = 0;
    f_vstr = FilterLFP(cfg_filter, CSC_vstr);
    
    %% Calculate LFP power
    cfg_power = [];
    IE = LFPpower(cfg_power, f_vstr);    % p_vstr is a TSD of power versus time of the filtered LFP trace.
    
    %% detect events
    cfg = [];
    cfg.method = 'percentile';
    cfg.threshold = cfg_out.percentile;
    cfg.dcn =  '>'; % return intervals where threshold is exceeded
    cfg.merge_thr = 0.05; % merge events closer than this
    cfg.minlen = 0.05; % minimum interval length
    
    lgp_evt = TSDtoIV(cfg,IE);
    
    

    event_IndicesStart = amplitude_peaks_index - numEventSamples/2; event_IndicesEnd = amplitude_peaks_index + numEventSamples/2;   % start and end indices either side of the peak value
    event_TimesStart = amplitude_peaks_times - cfg_out.eventDuration/2; event_TimesEnd = amplitude_peaks_times + cfg_out.eventDuration/2;    % start and end timestamps either side of the peak value
    
    %% Build IV
    eventStats = iv(event_TimesStart, event_TimesEnd);         % timestmaps for event start and end times
    eventStats.usr.istart = event_IndicesStart;                % tvec indices for start times
    eventStats.usr.iend = event_IndicesEnd;                    % tvec indices for end times
    eventStats.cfg.history.mfun = cat(1, eventStats.cfg.history.mfun, 'find_LFP_events.m'); % track what function worked on this data most recently
    eventStats.cfg.history.cfg = cat(1, eventStats.cfg.history.cfg, {cfg_out}); %  I don't think this is correct
    
    %% Arrange the LFP events into a data matrix
    disp('putting the data into array')
    tic
    nvars = 2; % change is to 3 if decide to use HIPP later
    nobs      = numEventSamples;   % number of observations per 'trial' (i.e. per gamma event)
    eventData = nan(nvars, nobs, length(amplitude_peaks_index));  % Channels x Datapoints in an event x Number of events
    disp(strcat('number of observations = ', num2str(nobs)))
    firstTimeStamps = amplitude_peaks_index - numEventSamples/2;
    firstGammaToUse = find(firstTimeStamps>0, 1, 'first');
    for iL = firstGammaToUse:length(amplitude_peaks_index);   % skip event 1 in case the gamma event 'starts' before the first timestamp.
        %     disp(iL);
        x = CSC_ofc.data(amplitude_peaks_index(iL)-(numEventSamples/2): amplitude_peaks_index(iL)+(numEventSamples/2));
        x = x(1:nobs);
        y = CSC_vstr.data(amplitude_peaks_index(iL)-(numEventSamples/2): amplitude_peaks_index(iL)+(numEventSamples/2));
        y = y(1:nobs);
        assert(size(x,2)==nobs);
        assert(size(y,2)==nobs);
        
        eventData(1,:,iL) = x;
        eventData(2,:,iL) = y;
    end
    toc
    
    %% plot it
    if cfg_out.doPlot == 1;
        clf
        hold on
        LoadExpKeys;
        plot(f_vstr.tvec, f_vstr.data);   % plot the filtered vStr LFP.
        plot(IE.tvec, IE.data, 'm')       % plot the amplitude envolope
        line([ExpKeys.TimeOnTrack ExpKeys.TimeOffTrack], [Amplitudecutoff Amplitudecutoff], 'linestyle', '--', 'Color', 'k', 'LineWidth', 1); % plot the threshold
        plot(CSC_vstr.tvec(amplitude_peaks_index), ones(1, length(amplitude_peaks_index)), 'k.', 'MarkerSize', cfg_out.MarkerSize)
        plot(amplitude_peaks_times, repmat(Amplitudecutoff, 1, length(amplitude_peaks_times)), 'c.', 'MarkerSize', cfg_out.MarkerSize)  % plot the event peaks
        c = axis;
        line([eventStats.tstart eventStats.tstart], [c(3) c(4)], 'color', 'g', 'LineWidth', 1);  % plot the event start times
        line([eventStats.tend eventStats.tend], [c(3) c(4)], 'color', 'r',  'LineWidth', 1);  % plot the event end times
        xlabel('Time (sec)', 'FontSize', 16)
        ylabel('Voltage (microvolts)', 'FontSize', 16)
        set(gca, 'FontSize', 16)
        title(ExpKeys.VSTRcsc);
    end
    %% save it
    if cfg_out.SaveIt == 1;
        SSN = GetSSN('SingleSession');
        fn = strcat(SSN, '-LFPevents-', eventName);
        save(fn, 'eventStats', 'eventData');
        disp(fn);
        disp('data saved');
    end
    popdir;
end

%% Criteria from mvdm lab papers
% must have 4 cycles
% must have 3 cycles > 50 microvolts amplitude
% variance score (variance/mean of cycle peaks and troughs)>1.5
% merge events that are separated by 50ms or less (same gamma type)
% Eric used a window of either 100ms or 3 cycles
% Julien used a window of 400ms

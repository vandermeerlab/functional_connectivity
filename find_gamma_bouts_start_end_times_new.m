function [iv_gammaTimes, X] = find_gamma_bouts_start_end_times_new(cfg_in)
% 2020-04-14. JJS. This function gets the start and end times of threshold crossing events for the filtered gamma LFP. 

% Inputs:
%   lowpass:            lowpass filter value (45 for g50, 70 for gamma50)
%   highpass:           highpass filter value (65 for g50, 90 for gamma80)
%   boutDuration:       how long you want the gamma bouts to be (in s). This parameter is chosen to make them all the same length, but in reality, gamma bouts
%                           will vary in duration.
% Outputs:
%   gammaIndicesStart:  tvec indices of the start of gamma events 
%   gammaIndicesEnd:    tvec indices of the end of gamma events
%   gammaTimesStart:    timestamps of the start of gamma events
%   gammaTimesEnd:      timestamps of the end of gamma events
%   X:                  nChannels x Samples x nGammaEvents           
                            % Channel 1 = OFC. Channel 2 = vstr     % can add Hipp later as Channel 3
                            % nSamples = number of data points in a gamma event (proportional to boutDuration)
                            % nGammaEvents = number of gamma events in the session 
                           
% Currently, function gets data (vStr and OFC) triggered from vStr gamma events only.                            

fd = FindFiles('*keys.m');
startSess = 1; endSess = length(fd);

cfg_def = [];
cfg_def.doPlot = 1;
cfg_def.MarkerSize = 20;
cfg_def.SaveIt = 0;
cfg_def.boutDuration = 0.1;  % Duration (in seconds) of the gamma events. This is somewhat arbitrary. Julien used 0.4 s, Eric used 0.1 s 
cfg_out = ProcessConfig(cfg_def,cfg_in);

for iSess = startSess:endSess;
    pushdir(fileparts(fd{iSess}));
    SSN = GetSSN('SingleSession'); disp(SSN);
    
    cfg_prep = [];
    cfg_prep.decimateByFactor = 2; % Downsampling the data. Default here is 2 (from 2kHz to 1kHz).
    cfg_prep.detrend = 1;   % for removing slow DC shifts in voltage
    cfg_prep.diffdata = 0;  % for removing autocorrelation in the time series
    cfg_prep.doRestrict = 1; % Restirct CSC to track time only.
    cfg_prep.hippflag = 0;  % Do we want to Load a hippocampal CSC? (for sessions with electrode in fissure). 0 = no. 1 = yes.
    [ofc, vstr, ~, ~] = prepCSCs_new(cfg_prep);   % can change this later to return hippocampus LFP also
    fs = 1/median(diff(vstr.tvec));  % sampling rate of the LFP
    
    %% Filter in low gamma band
    cfg_filter = [];
    cfg_filter.f = [45 65];  % frequency band to filter in. So far, just using low gamma events [45 65]
    cfg_filter.display_filter = 0;
    f_vstr = FilterLFP(cfg_filter,vstr);
    
    %% Calculate gamma power
    cfg_power = [];
    IE = LFPpower(cfg_power,f_vstr);    % p_vstr is a TSD of power versus time of the filtered LFP trace.
    percentile = .95;                   % Find peaks in the amplitude envelope that are above the 95% level. 
    IEsort = sort(IE.data);
    percentile_index = percentile*length(IEsort);
    Amplitudecutoff = IEsort(round(percentile_index));
    TimeStepsToChopOFF = 2*round(cfg_out.boutDuration*fs); % How many samples to ignore at the end of the session...so that the last gamma event doesn't exceed the last timestep
    numBoutSamples = round(cfg_out.boutDuration*fs/2);  % number of samples on either side of the peak of a gamma event to include 
    if mod(numBoutSamples,2) == 1;
        numBoutSamples = numBoutSamples-1;
    end
    
    %% Find Threshold Crossings 
    % Get timestamps for the first instance of each threshold crossing for the gamma amplitude envelope
    disp('finding threshold crossings')
    [~, amplitdue_index] = find(IE.data>=Amplitudecutoff);
    event_diff = cat(2, 1, diff(amplitdue_index));
    transition_points = event_diff>10;
    event_boundaries_index = amplitdue_index(transition_points==1);
    % event_boundaries = IE.T(event_boundaries_index);
    
    %% By event INDICES
    % demarcate the gamma events in a matrix.
    lasttouse = length(IE.data)-TimeStepsToChopOFF;
    event_boundaries_touse = event_boundaries_index<lasttouse;
    event_boundaries_index = event_boundaries_index(event_boundaries_touse); % to remove cases where gamma event occurs in last 100ms of the recording
    for iEvent = 1:length(event_boundaries_index);
        t0 = event_boundaries_index(iEvent);
        datatouse = IE.data(t0: t0 + 100); % look in the next 100ms.   * this is approximate. 100 samples at ~ 1 sample/millisecond
        [~, c] = max(datatouse);
        amplitude_peaks_index(iEvent) = t0 + c;
        amplitude_peaks_times(iEvent) = IE.tvec(t0 + c);
    end
    amplitude_peaks_index = amplitude_peaks_index';
    gammaIndicesStart = amplitude_peaks_index - numBoutSamples; gammaIndicesEnd = amplitude_peaks_index + numBoutSamples;
    gammaTimesStart = amplitude_peaks_times - numBoutSamples; gammaTimesEnd = amplitude_peaks_times + numBoutSamples;
    
    %% Build IV 
    iv_gammaTimes = iv(gammaTimesStart, gammaTimesEnd); 
  
    
    
    %% Pull out the data
    disp('putting the data into array')
    tic
    nvars = 2; % change is to 3 if decide to use HIPP later
    nobs      = numBoutSamples*2;   % number of observations per 'trial' (i.e. per gamma event)
    X = nan(nvars, nobs, length(amplitude_peaks_index));  % for nvars, ofc = 1, vstr = 2;
    disp(strcat('number of observations = ', num2str(nobs)))
    firstTimeStamps = amplitude_peaks_index - numBoutSamples;
    firstGammaToUse = find(firstTimeStamps>0, 1, 'first');
    for iL = firstGammaToUse:length(amplitude_peaks_index);   % skip event 1 in case the gamma event 'starts' before the first timestamp.
        %     disp(iL);
        x = ofc.data(amplitude_peaks_index(iL)-(numBoutSamples): amplitude_peaks_index(iL)+(numBoutSamples));
        x = x(1:nobs);
        y = vstr.data(amplitude_peaks_index(iL)-(numBoutSamples): amplitude_peaks_index(iL)+(numBoutSamples));
        y = y(1:nobs);
        assert(size(x,2)==nobs);
        assert(size(y,2)==nobs);
        
        X(1,:,iL) = x;
        X(2,:,iL) = y;
    end
    toc
    
    %% plot it
    if cfg_out.doPlot == 1;
        clf
        hold on
        LoadExpKeys;
        plot(f_vstr.tvec, f_vstr.data);   % plot the filtered vStr LFP.
        plot(IE.tvec, IE.data, 'r')       % plot the amplitude envolope 
        line([ExpKeys.TimeOnTrack ExpKeys.TimeOffTrack], [Amplitudecutoff Amplitudecutoff], 'linestyle', '--', 'Color', 'k', 'LineWidth', 1);
        plot(vstr.tvec(amplitude_peaks_index), ones(1,length(amplitude_peaks_index)), 'k.', 'MarkerSize', cfg_out.MarkerSize)
        plot(amplitude_peaks_times, repmat(Amplitudecutoff,1,length(amplitude_peaks_times)), 'g.', 'MarkerSize', cfg_out.MarkerSize)
        xlabel('Time (sec)', 'FontSize', 16)
        ylabel('Voltage (microvolts)', 'FontSize', 16)
        set(gca, 'FontSize', 16)
        title(ExpKeys.VSTRcsc);
    end
    %% save it
    if cfg_out.SaveIt == 1;
        if lowpass == 45; gammatype = 'low'; end
        if lowpass == 70; gammatype = 'high'; end
        SSN = GetSSN('SingleSession');
        fn = strcat(SSN, '-gamma_bouts-', gammatype);
        save(fn, 'gammaIndicesStart', 'gammaIndicesEnd', 'gammaTimesStart', 'gammaTimesEnd', 'X', 'params');
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

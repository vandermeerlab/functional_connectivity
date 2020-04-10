function [gammaIndicesStart, gammaIndicesEnd, gammaTimesStart, gammaTimesEnd, X] = find_gamma_bouts_start_end_times_new(cfg_in)
% Inputs:
%   lowpass:        lowpass filter value (45 for g50, 70 for g80)
%   highpass:       highpass filter value (65 for g50, 90 for g80)
%   boutDuration:   how long you want the gamma bouts to be (in s). This parameter is chosen to make them all the same length, but in reality, gamma bouts
%                       will vary in duration.
% Outputs:
%   amplitude_peaks = timestamps of gamma amplitude peaks
%
fd = FindFiles('*keys.m');
startSess = 1; endSess = length(fd);

cfg_def.MarkerSize = 20;
cfg_def.SaveIt = 1;
cfg_def.doPlot = 0;
cfg_def.decimateFactor = 2;
cfg_def.doRestrict = 1;
cfg_def.FiltOrder = 256;
cfg_def.lowpass = 45;
cfg_def.highpass = 65;
cfg_def.boutDuration = 0.1;
master_cfg = ProcessConfig(cfg_def,cfg_in);

for iSess = startSess:endSess;
    pushdir(fileparts(fd{iSess}));
    SSN = GetSSN('SingleSession'); disp(SSN);
    
    [ofc, vstr, ~, ~] = prepCSCs_new(master_cfg);   % can change this later to return hippocampus LFP also 
    
    %% NEW
    % filter in low gamma band
    cfg = [];
    cfg.f = [master_cfg.lowpass master_cfg.highpass];
    cfg.display_filter = 0;   
    fvstr = FilterLFP(cfg,vstr);
    
        %%  OLD
%     disp('calculating power')
%     % [IF, IA, IP, CSC0, IE]=InstSig(vstr,45,65,256);
%     [~, ~, ~, CSC0, IE] = InstSig(vstr,lowpass,highpass,FiltOrder);
%     percentile = .95;
%     IEsort = sort(IE.data);
%     percentile_index = percentile*length(IEsort);
%     Amplitudecutoff = IEsort(round(percentile_index));
%     TimeStepsToChopOFF = 4*round(boutDuration*fs);
%     numBoutSamples = round(boutDuration*fs/2);  % number of samples on either side of the peak of a gamma event
%     if mod(numBoutSamples,2) == 1;
%         numBoutSamples = numBoutSamples-1;
%     end  
    %%
    disp('finding threshold crossings')
    % Get timestamps for the first instance of each threshold crossing for the gamma amplitude envelope
    [amplitdue_index,~] = find(IE.data>=Amplitudecutoff);
    event_diff = cat(1, 1, diff(amplitdue_index));
    transition_points = event_diff>10;
    event_boundaries_index = amplitdue_index(transition_points==1);
    % event_boundaries = IE.T(event_boundaries_index);
    
    disp('finding peak of candidate events')
    %% By event INDICES
    % demarcate the gamma events in a matrix.
    lasttouse = length(IE.D)-TimeStepsToChopOFF;
    event_boundaries_touse = event_boundaries_index<lasttouse;
    event_boundaries_index = event_boundaries_index(event_boundaries_touse); % to remove cases where gamma event occurs in last 100ms of the recording
    for iEvent = 1:length(event_boundaries_index);
        t0 = event_boundaries_index(iEvent);
        datatouse = IE.D(t0: t0 + 100); % look in the next 100ms.   * this is approximate. 100 samples at ~ 1 sample/millisecond
        [~, c] = max(datatouse);
        amplitude_peaks_index(iEvent) = t0 + c;
        amplitude_peaks_times(iEvent) = IE.T(t0 + c);
    end
    amplitude_peaks_index = amplitude_peaks_index';
    gammaIndicesStart = amplitude_peaks_index - numBoutSamples; gammaIndicesEnd = amplitude_peaks_index + numBoutSamples;
    gammaTimesStart = amplitude_peaks_times - numBoutSamples; gammaTimesEnd = amplitude_peaks_times + numBoutSamples;
    
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
        x = ofc.D(amplitude_peaks_index(iL)-(numBoutSamples): amplitude_peaks_index(iL)+(numBoutSamples));
        x = x(1:nobs);
        y = vstr.D(amplitude_peaks_index(iL)-(numBoutSamples): amplitude_peaks_index(iL)+(numBoutSamples));
        y = y(1:nobs);
        assert(size(x,1)==nobs);
        assert(size(y,1)==nobs);
        
        X(1,:,iL) = x;
        X(2,:,iL) = y;
    end
    toc
    
    %% plot it
    if doPlot == 1;
        clf
        hold on
        EvalKeys;
        plot(CSC0.range, CSC0.data);
        plot(IE.range, IE.data, 'r')
        line([ExpKeys.TimeOnTrack ExpKeys.TimeOffTrack], [Amplitudecutoff Amplitudecutoff], 'linestyle', '--', 'Color', 'k', 'LineWidth', 1);
        plot(vstr.T(amplitude_peaks_index), ones(1,length(amplitude_peaks_index)), 'k.', 'MarkerSize', MarkerSize)
        plot(amplitude_peaks_times, repmat(Amplitudecutoff,1,length(amplitude_peaks_times)), 'g.', 'MarkerSize', MarkerSize)
    end
    %% save it
    if SaveIt == 1;
        if lowpass == 45; gammatype = 'low'; end
        if lowpass == 70; gammatype = 'high'; end
        SSN = GetSSN('SingleSession');
        fn = strcat(SSN, '-gamma_bouts-', gammatype);
        save(fn, 'gammaIndicesStart', 'gammaIndicesEnd', 'gammaTimesStart', 'gammaTimesEnd', 'X', 'params');
        disp(fn);
        disp('data saved');
    end
    
    %     clear gammaIndicesStart; clear gammaIndicesEnd; clear gammaTimesStart; clear gammaTimesEnd; clear params; clear amplitude_peaks_index; clear amplitude_peaks_times;
    popdir;
end


%% Criteria from mvdm lab papers
% must have 4 cycles
% must have 3 cycles > 50 microvolts amplitude
% variance score (variance/mean of cycle peaks and troughs)>1.5
% merge events that are separated by 50ms or less (same gamma type)
% Eric used a window of either 100ms or 3 cycles
% Julien used a window of 400ms

function [eventStats, eventData] = find_gamma_bouts_start_end_times_new_single_session(cfg_in, eventName)
% 2020-04-14. JJS. This function gets the start and end times of threshold crossing events for the filtered gamma LFP.
% Operates on a single session and a single frequency band. 

% Inputs:
%   cfg_in:                         config input
%   eventName                       name (string) that the file will be saved as; for instance 'gamma50'

% Outputs:
%   eventStats:                     IV structure with event start and end times

%   eventData:                      nChannels x Samples x nEvents
%                                       Typically, channel 1 = OFC. channel 2 = vStr     % can add Hipp later as Channel 3
%                                       nSamples = number of data points in an LFP event (proportional to eventDuration)
%                                       nEvents = number of events in the session

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
percentile = cfg_out.percentile;                   % Find peaks in the amplitude envelope that are above the 95% level.
IEsort = sort(IE.data);
percentile_index = percentile*length(IEsort);
Amplitudecutoff = IEsort(round(percentile_index));
TimeStepsToChopOFF = round(cfg_out.eventDuration*fs)+1; % How many samples to ignore at the end of the session...so that the last gamma event doesn't exceed the last timestep
numEventSamples = round(cfg_out.eventDuration*fs);          % number of samples in an event. Divide by 2 to get numSamples on either side of peak.
if mod(numEventSamples/2,2) == 1;
    numEventSamples = numEventSamples-1;
end
%% Find Threshold Crossings
% Get timestamps for the first instance of each threshold crossing for the filtered LFP envelope
disp('finding threshold crossings')
[~, amplitdue_index] = find(IE.data>=Amplitudecutoff);          % get the set of data points where the LFP is greater than the percentile cuttoff
event_diff = cat(2, 1, diff(amplitdue_index));                  % diff(amplitude_index) returns ones for adjacent points (part of a contiguous event), and larger values for the spaces between events
transition_points = event_diff>cfg_out.minIEI*fs;               % find where the event transitions occur. This is a logical with ones equal to the transition points. Must be farther apart than length of an event (eventDuration). Otherwise, will be skipped.
event_boundaries_index = amplitdue_index(transition_points==1); % index of times where events BEGIN

%% By event INDICES
% demarcate the LFP events in a matrix.
lasttouse = length(IE.data)-TimeStepsToChopOFF;                 % find the last time point to consider or inclusion
event_boundaries_touse = event_boundaries_index<lasttouse;      % use all candidate events but those falling past teh cutoff
event_boundaries_index = event_boundaries_index(event_boundaries_touse); % to remove cases where LFP event occurs in last bit of the recording < eventDuration
for iEvent = 1:length(event_boundaries_index);
    t0 = event_boundaries_index(iEvent);
    datatouse = IE.data(t0: t0 + numEventSamples);              % time indices that consititue the LFP event
    [~, c] = max(datatouse);                                    % find which value is the peak of the LFP event
    amplitude_peaks_index(iEvent) = t0 + c;                     % get the index
    amplitude_peaks_times(iEvent) = IE.tvec(t0 + c);            % get the timestamp
end
amplitude_peaks_index = amplitude_peaks_index';
event_IndicesStart = amplitude_peaks_index - numEventSamples/2; event_IndicesEnd = amplitude_peaks_index + numEventSamples/2;   % start and end indices either side of the peak value
event_TimesStart = amplitude_peaks_times - cfg_out.eventDuration/2; event_TimesEnd = amplitude_peaks_times + cfg_out.eventDuration/2;    % start and end timestamps either side of the peak value

%% Build IV
eventStats = iv(event_TimesStart, event_TimesEnd);         % timestmaps for event start and end times
eventStats.usr.istart = event_IndicesStart;                % tvec indices for start times
eventStats.usr.iend = event_IndicesEnd;                    % tvec indices for end times
eventStats.cfg.history.mfun = cat(1, eventStats.cfg.history.mfun, 'find_LFP_events.m'); % track what function worked on this data most recently
eventStats.cfg.history.cfg = cfg_out; %  I don't think this is correct
%     eventStats.cfg.history.cfg = cat(1, eventStats.cfg.history.cfg, {cfg_out}); %  I don't really understand this line
eventStats.usr.fs = fs; % make a note of the sampling rate

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
    x = CSC_ofc.data(amplitude_peaks_index(iL)-(numEventSamples/2): amplitude_peaks_index(iL)+(numEventSamples/2));   % OFC
    x = x(1:nobs);
    y = CSC_vstr.data(amplitude_peaks_index(iL)-(numEventSamples/2): amplitude_peaks_index(iL)+(numEventSamples/2));  % vStr
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
%     pause;
end
%% save it
if cfg_out.SaveIt == 1;
    SSN = GetSSN('SingleSession');
    fn = strcat(SSN, '-LFPevents-', eventName);
    save(fn, 'eventStatsToSave', 'eventDataToSave');
    disp(fn);
    disp('data saved');
end


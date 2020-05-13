function [X] = assign_LFP_events(eventName, varargin)
%04-2020. JJS. Assign LFP events as occuring in one of several different task epochs.

% Inputs:                       X - structure with fields that include
%                                   X.eventStats.usr.istart   - tvec start Indices for LFP events
%                                   X.eventStats.usr.iend     - tvec end Indices for LFP events
%
% Outputs:
%
X = [];
% skipSess = NaN;
doPlot = 1;
doSave = 1;
MarkerSize = 15;
rewardStart = 1;
rewardEnd = 3;
process_varargin(varargin);

if doPlot == 1;
    plotPosition = 1;
else
    plotPosition = 0;
end

LoadExpKeys;
sd = DDinit_new;

SSN = GetSSN('SingleSession');
fn = strcat(SSN, '-XCorr-', eventName, '.mat');

if exist(fn, 'file') == 2;           % Load LFP events. If file is not present, skip session. 
    load(fn, 'X');
       
    %% Load the Position data
    if exist('VT1.Nvt', 'file') == 2;
        pos_tsd = LoadPos([]);
        % elseif exist(strcat(SSN, '-vt.mat'), 'file');
        %     fn = strcat(SSN, '-vt.mat');
        %     load(fn);
    end
    
    tstart = X.eventStats.tstart; %#ok<NODEF>
    
    %% Find all gamma events during the task 
    after =  ExpKeys.TimeOnTrack < tstart;
    before = tstart < ExpKeys.TimeOffTrack;  
    TaskEvents = before & after & X.keep' == 1;
    X.TaskEvents = TaskEvents';
    if plotPosition; mazeTask_idx = TSD_getidx(pos_tsd, X.eventStats.tstart(TaskEvents), X.eventStats.tend(TaskEvents)); end
    
    %% Find all gamma events on pedastal 
    PRE =  ExpKeys.TimeOnTrack > tstart;
    POST = tstart < ExpKeys.TimeOffTrack;  
    Pedastal = PRE & POST & X.keep' == 1;
    X.Pedastal = Pedastal';
    
    %% Find gamma events during PRE-run pedastal period
    PreRun_LFPevents = tstart < ExpKeys.TimeOnTrack  & X.keep' == 1;   % these are indices for LFP events in X.eventStats that fall within the reward window. These are in terms of event order, NOT cheetah timestamps.
    
    PreRun_LFPevents = PreRun_LFPevents(~isnan(PreRun_LFPevents));  % there are rare NANs for sd.FeederFire times. Cannot index NaNs below
    X.PreRun_LFPevents = PreRun_LFPevents';
    if plotPosition; mazePreRun_idx = TSD_getidx(pos_tsd, X.eventStats.tstart(PreRun_LFPevents), X.eventStats.tend(PreRun_LFPevents)); end
    
    %% Find gamma events during POST-run pedastal period
    PostRun_LFPevents = tstart > ExpKeys.TimeOffTrack  & X.keep' == 1;   % these are indices for LFP events in X.eventStats that fall within the reward window. These are in terms of event order, NOT cheetah timestamps.
    
    PostRun_LFPevents = PostRun_LFPevents(~isnan(PostRun_LFPevents));  % there are rare NANs for sd.FeederFire times. Cannot index NaNs below
    X.PostRun_LFPevents = PostRun_LFPevents';
    if plotPosition; mazePostRun_idx = TSD_getidx(pos_tsd, X.eventStats.tstart(PostRun_LFPevents), X.eventStats.tend(PostRun_LFPevents)); end
    
    %% Find gamma events proximal to reward delivery
    Reward_LFPevents = false(length(X.keep),1);
    for iLap = 1:length(sd.FeederTimes);
        idx = tstart > sd.FeederTimes(iLap) + rewardStart & tstart < sd.FeederTimes(iLap) + rewardEnd & X.keep' == 1;
        Reward_LFPevents = idx | Reward_LFPevents;
    end
    X.Reward_LFPevents = Reward_LFPevents';
    if plotPosition; mazeReward_idx = TSD_getidx(pos_tsd, X.eventStats.tstart(Reward_LFPevents), X.eventStats.tend(Reward_LFPevents)); end
    
    %% Find gamma events that occur during VTE choice point passes
    VTE_LFPevents = false(length(X.keep),1);
    for iLap = 1:length(sd.VTElaps);
        idx = tstart > sd.EnteringCPTime(sd.VTElaps(iLap)) & tstart < sd.ExitingCPTime(sd.VTElaps(iLap)) & X.keep' == 1;
        VTE_LFPevents = idx | VTE_LFPevents;  % these are indices for LFP events in X.eventStats that fall within the reward window. These are in terms of event order, NOT cheetah timestamps.
    end
    X.VTE_LFPevents = VTE_LFPevents';
    if plotPosition; mazeVTE_idx = TSD_getidx(pos_tsd, X.eventStats.tstart(VTE_LFPevents), X.eventStats.tend(VTE_LFPevents)); end
    
    %% Find gamma events that occur during non-VTE choice point passes
    nonVTE_LFPevents = false(length(X.keep),1);
    for iLap = 1:length(sd.nonVTElaps);
        idx = tstart > sd.EnteringCPTime(sd.nonVTElaps(iLap)) & tstart < sd.ExitingCPTime(sd.nonVTElaps(iLap)) & X.keep' == 1;
        nonVTE_LFPevents = idx | nonVTE_LFPevents;  % these are indices for LFP events in X.eventStats that fall within the reward window. These are in terms of event order, NOT cheetah timestamps.
    end
    X.nonVTE_LFPevents = nonVTE_LFPevents';
    if plotPosition; mazenonVTE_idx = TSD_getidx(pos_tsd, X.eventStats.tstart(nonVTE_LFPevents), X.eventStats.tend(nonVTE_LFPevents)); end
    
    %% Find gamma events that occur during Feeder Approach
    Approach_LFPevents = false(length(X.keep),1);
    for iLap = 1:length(sd.EnteringZoneTime);
        idx = tstart > sd.EnteringZoneTime(iLap) & tstart < sd.EnteringZoneTime(iLap) + 1 & X.keep' == 1;  % 0.75sec gets you close to but not yet at the feeder
        Approach_LFPevents = idx | Approach_LFPevents;  % these are indices for LFP events in X.eventStats that fall within the reward window. These are in terms of event order, NOT cheetah timestamps.
    end
    X.Approach_LFPevents = Approach_LFPevents';
    if plotPosition; mazeApproach_idx = TSD_getidx(pos_tsd, X.eventStats.tstart(Approach_LFPevents), X.eventStats.tend(Approach_LFPevents)); end
    
    %% Find gamma events that occur during delay side Reward Waiting
    Waiting_LFPevents = false(length(X.keep),1);
    for iLap = 1:length(sd.D);
        idx = tstart > sd.D(iLap) & tstart -2 < sd.D(iLap) & X.keep' == 1;  % 0.75sec gets you close to but not yet at the feeder
        Waiting_LFPevents = idx | Waiting_LFPevents;  % these are indices for LFP events in X.eventStats that fall within the reward window. These are in terms of event order, NOT cheetah timestamps.
    end
    X.Waiting_LFPevents = Waiting_LFPevents';
    if plotPosition; mazeWaiting_idx = TSD_getidx(pos_tsd, X.eventStats.tstart(Waiting_LFPevents), X.eventStats.tend(Waiting_LFPevents)); end
    
    % %% Find gamma events that occur start maze segment
    % SoM_LFPevents = NaN;
    % for iLap = 1:length(sd.D);
    %     idx = find(tstart > sd.EnteringSoMTime(iLap) & tstart < sd.ExitingSoMTime(iLap) & X.keep' == 1);
    %     SoM_LFPevents = vertcat(SoM_LFPevents, idx);  % these are indices for LFP events in X.eventStats that fall within the reward window. These are in terms of event order, NOT cheetah timestamps.
    % end
    % SoM_LFPevents = SoM_LFPevents(~isnan(SoM_LFPevents));
    % X.SoM_LFPevents = SoM_LFPevents;
    % % X.SoM_LFPindices = X.eventStats.usr.istart(SoM_LFPevents);
    % % X.SoM_LFPtimes = X.eventStats.tstart(SoM_LFPevents);
    % if plotPosition; mazeSoM_idx = TSD_getidx(pos_tsd, X.eventStats.tstart(SoM_LFPevents), X.eventStats.tend(SoM_LFPevents)); end
    
    disp(strcat('fract covered events =', num2str((sum(X.Reward_LFPevents) + sum(X.VTE_LFPevents) + ...
        sum(X.nonVTE_LFPevents) + sum(X.Approach_LFPevents) + sum(X.Waiting_LFPevents) ...
        + sum(X.PreRun_LFPevents) + sum(X.PostRun_LFPevents))/sum(X.keep))));
    
    %% Plot the results to see if it looks right
    if doPlot == 1;
        figure(100)
        clf
        subplot(1,2,1) 
        plot(-pos_tsd.data(2,:), -pos_tsd.data(1,:), '.', 'MarkerSize', 1, 'color', [.5 .5 .5]); hold on
        plot(-pos_tsd.data(2,mazeTask_idx), -pos_tsd.data(1,mazeTask_idx), 'b.', 'MarkerSize', MarkerSize)
        legend('all position samples', 'all task events', 'Location', 'SouthEast')
        title(strcat(SSN, '-', eventName), 'FontSize', 15)
        
        subplot(1,2,2) 
        plot(-pos_tsd.data(2,:), -pos_tsd.data(1,:), '.', 'MarkerSize', 1, 'color', [.5 .5 .5]); hold on
        plot(-pos_tsd.data(2,mazePreRun_idx), -pos_tsd.data(1,mazePreRun_idx), 'k.', 'MarkerSize', MarkerSize) % pre-run, pedastal
        %     plot(-pos_tsd.data(2,mazeSoM_idx), -pos_tsd.data(1,mazeSoM_idx), 'c.', 'MarkerSize', MarkerSize) % Start of Maze
        plot(-pos_tsd.data(2,mazenonVTE_idx), -pos_tsd.data(1,mazenonVTE_idx), 'y.', 'MarkerSize', MarkerSize)  % nonVTE passes
        plot(-pos_tsd.data(2,mazeVTE_idx), -pos_tsd.data(1,mazeVTE_idx), 'b.', 'MarkerSize', MarkerSize)   % VTE passes
        plot(-pos_tsd.data(2,mazeApproach_idx), -pos_tsd.data(1,mazeApproach_idx), 'r.', 'MarkerSize', MarkerSize) % Feeder Approach
        plot(-pos_tsd.data(2,mazeWaiting_idx), -pos_tsd.data(1,mazeWaiting_idx), 'm.', 'MarkerSize', MarkerSize)  % Delay side Waiting
        plot(-pos_tsd.data(2,mazeReward_idx), -pos_tsd.data(1,mazeReward_idx), 'g.', 'MarkerSize', MarkerSize)  % Reward site LFP events
        plot(-pos_tsd.data(2,mazePostRun_idx), -pos_tsd.data(1,mazePostRun_idx), '.', 'MarkerSize', MarkerSize, 'color', [.5 1 1]) % post-run, pedastal
        %     pause
        legend('all position samples', 'pre-run', 'nonVTE passes', 'VTE passes', 'approach', 'waiting', 'reward consumption', 'post-run', 'Location', 'SouthEast')
        title(strcat(SSN, '-', eventName), 'FontSize', 15)
    end
    
    if doSave == 1;
        assert(~isempty(X))
        fn = strcat(SSN, '-XCorr-', eventName);
        save(fn, 'X');
        disp('LFP epochs saved')
    end
       
else
    disp('LFP events not found. Skipping session')
%     [~, skipSess, ~] = fileparts(pwd);
end









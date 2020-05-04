function [X] = assign_LFP_events(eventName, varargin)
%04-2020. JJS. Assign LFP events as occuring in one of several different task epochs.

% Inputs:                       X - structure with fields that include
%                                   X.eventStats.usr.istart   - tvec start Indices for LFP events
%                                   X.eventStats.usr.iend     - tvec end Indices for LFP events
%
% Outputs:                          X.usr.
%
process_varargin(varargin);

LoadExpKeys;
sd = DDinit_new;

SSN = GetSSN('SingleSession');
fn = strcat(SSN, '-XCorr-', eventName);
load(fn, 'X');

% LFPindicesStartToUse = X.eventStats.usr.istart(X.keep ==1);
% LFPindicesEndToUse = X.eventStats.usr.istart(X.keep ==1);

tstart = X.eventStats.tstart(X.keep ==1); %#ok<NODEF>
% LFPtimesEndToUse = X.eventStats.tend(X.keep ==1);

%% Find gamma events proximal to reward and those proximal to ChoicePoint passes
Reward_LFPindex = NaN;
for iLap = 1:length(sd.FeederTimes);
    indicesToUse = find(tstart > sd.FeederTimes(iLap) & tstart < sd.FeederTimes(iLap)+1);  % three seconds after FeederFire
    Reward_LFPindex = vertcat(Reward_LFPindex, indicesToUse);
end
Reward_LFPindex = Reward_LFPindex(~isnan(Reward_LFPindex));
X.RewardLFPindex = Reward_LFPindex;  
X.RewardLFPtimes = X.eventStats.tstart(Reward_LFPindex);

CP_LFPindex = NaN;
for iLap = 1:length(sd.EnteringCPTime);
    indicesToUse = find(tstart > sd.EnteringCPTime(iLap) & tstart < sd.ExitingCPTime(iLap));
    CP_LFPindex = vertcat(CP_LFPindex, indicesToUse);
end
CP_LFPindex = CP_LFPindex(~isnan(CP_LFPindex));
X.CP_LFPindex = CP_LFPindex; 
X.CP_LFPtimes = X.eventStats.tstart(CP_LFPindex);

end


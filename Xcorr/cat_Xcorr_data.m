function [Xcat] = cat_Xcorr_data(eventName, varargin)
%2020-05-09. JJS. Concatenate amplitude Xcorr data that is saved within each session folder.
%   Gathers Xcorr values from each session that is saved in the structure 'X' by amp_crosscorr_on_gamma_events_new.m

fd = FindFiles('*keys.m');
startSess = 1;
endSess = length(fd);
process_varargin(varargin);

Xcat.crosscorr = [];
Xcat.max_crosscorr_lag = [];
Xcat.keep = [];
% Xcat.Allpedastal = [];
Xcat.AllTask = [];
Xcat.PreRun_LFPevents = [];
Xcat.PostRun_LFPevents = [];
Xcat.Reward_LFPevents = [];
Xcat.VTE_LFPevents = [];
Xcat.nonVTE_LFPevents = [];
Xcat.Approach_LFPevents = [];
Xcat.Waiting_LFPevents = [];

for iSess = startSess:endSess;
    pushdir(fileparts(fd{iSess}));
    SSN = GetSSN('SingleSession');
    disp(SSN);
    fn = strcat(SSN, '-XCorr-', eventName);
    
    if exist(strcat(fn, '.mat'), 'file') == 2;
        load(fn, 'X');
        
        Xcat.lags = X.lags;
        Xcat.crosscorr = cat(2, Xcat.crosscorr, X.crosscorr);
        Xcat.max_crosscorr_lag = cat(2, Xcat.max_crosscorr_lag, X.max_crosscorr_lag);
        Xcat.keep = cat(2, Xcat.keep, X.keep);
        
        %     Xcat.Allpedastal = cat(2, Xcat.Allpedastal, X.Pedastal);  % note: all of the epochs from X already have X.keep incorporated
        Xcat.AllTask = cat(2, Xcat.AllTask, X.TaskEvents);
        Xcat.PreRun_LFPevents = cat(2, Xcat.PreRun_LFPevents, X.PreRun_LFPevents);
        Xcat.PostRun_LFPevents = cat(2, Xcat.PostRun_LFPevents, X.PostRun_LFPevents);
        Xcat.Reward_LFPevents = cat(2, Xcat.Reward_LFPevents, X.Reward_LFPevents);
        Xcat.VTE_LFPevents = cat(2, Xcat.VTE_LFPevents, X.VTE_LFPevents);
        Xcat.nonVTE_LFPevents = cat(2, Xcat.nonVTE_LFPevents, X.nonVTE_LFPevents);
        Xcat.Approach_LFPevents = cat(2, Xcat.Approach_LFPevents, X.Approach_LFPevents);
        Xcat.Waiting_LFPevents = cat(2, Xcat.Waiting_LFPevents, X.Waiting_LFPevents);
        
        popdir;
    end
end

Xcat.keep = logical(Xcat.keep);
% Xcat.Allpedastal = [];
Xcat.AllTask = logical(Xcat.AllTask);
Xcat.PreRun_LFPevents = logical(Xcat.PreRun_LFPevents);
Xcat.PostRun_LFPevents = logical(Xcat.PostRun_LFPevents);
Xcat.Reward_LFPevents = logical(Xcat.Reward_LFPevents);
Xcat.VTE_LFPevents = logical(Xcat.VTE_LFPevents);
Xcat.nonVTE_LFPevents = logical(Xcat.nonVTE_LFPevents);
Xcat.Approach_LFPevents = logical(Xcat.Approach_LFPevents);
Xcat.Waiting_LFPevents = logical(Xcat.Waiting_LFPevents);


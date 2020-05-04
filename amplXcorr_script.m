function amplXcorr_script(varargin)
% Runs through the functions to get LFP event times and calculate XCorr for a single session.

fd = FindFiles('*keys.m');
startSess = 1; 
endSess = length(fd);
process_varargin(varargin);

for iSess = startSess:endSess;
    pushdir(fileparts(fd{iSess}));
    SSN = GetSSN('SingleSession');
    disp(SSN);
    
    [CSC_ofc, CSC_vstr, ~, ~, cfg_csc] = prepCSCs_new([]);
    
    [eventStats, eventData] = find_gamma_bouts_start_end_times_new_single_session([], CSC_ofc, CSC_vstr, 'highgamma');
    
    [X, cfg_amp] = amp_crosscorr_on_gamma_events_new([], CSC_ofc, CSC_vstr, eventStats, eventData, 'highgamma', 'doPlotKeeps', 0);
    
end

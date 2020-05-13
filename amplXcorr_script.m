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
    
    [eventStats, eventData] = find_LFP_events([], CSC_ofc, CSC_vstr, 'beta');
    
    [X, cfg_amp] = Xcorr_on_LFP_events([], CSC_ofc, CSC_vstr, eventStats, eventData, 'beta', 'doPlotKeeps', 0);
    
end

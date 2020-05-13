function assign_LFP_events_all_sessions(eventName, varargin)

fd = FindFiles('*keys.m');
startSess = 1;
endSess = length(fd);
process_varargin(varargin);

for iSess = startSess:endSess;
    pushdir(fileparts(fd{iSess}));
    SSN = GetSSN('SingleSession');
    disp(SSN);
    [~] = assign_LFP_events(eventName, 'doPlot', 0, 'doSave', 1);
    popdir;
end

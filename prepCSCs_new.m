function [ofc, vstr, hipp, fs, cfg_csc] = prepCSCs_new(cfg_in)

%2014-09-10. JJS. Loads and prepares the ofc and vstr CSC files that are
%specified as 'good' in the keys file for analysis of granger causality for
%a single session.
%2020-04-06. JJS. Modifications to work w/ MvdMlab codeset.
cfg_def = [];
cfg_def.VoltageConvFactor = 10^6;
cfg_def.decimateByFactor = 2; % Downsampling the data. Default here is 2 (from 2kHz to 1kHz).
cfg_def.detrend = 1;   % for removing slow DC shifts in voltage
cfg_def.diffdata = 0;  % for removing autocorrelation in the time series
cfg_def.doRestrict = 1; % Restirct CSC to track time only.
cfg_def.hippflag = 0;  % Do we want to Load a hippocampal CSC? (for sessions with electrode in fissure). 0 = no. 1 = yes.
cfg_csc = ProcessConfig(cfg_def,cfg_in);

LoadExpKeys;

%% Load the CSCs
cfg_csc.fc = {ExpKeys.OFCcsc};
ofc = LoadCSC(cfg_csc);
ofc_dt = median(diff(ofc.tvec));

cfg_csc.fc = {ExpKeys.VSTRcsc};
vstr = LoadCSC(cfg_csc);
vstr_dt = median(diff(vstr.tvec));

fs = 1/median(diff(ofc.tvec)); 

if strcmp('NaN', ExpKeys.HIPPcsc) == 0 && cfg_csc.hippflag == 1;
    cfg_csc.fc = ExpKeys.HIPPcsc;
    [hipp] = LoadCSC(cfg_csc);
    hipp_dt = median(diff(hipp.tvec));
else
    hipp = [];
end


%% restrict to time on track
if cfg_csc.doRestrict == 1;
    ofc = restrict(ofc, ExpKeys.TimeOnTrack, ExpKeys.TimeOffTrack);
    vstr = restrict(vstr, ExpKeys.TimeOnTrack, ExpKeys.TimeOffTrack);
end

if cfg_csc.doRestrict == 1 && cfg_def.hippflag == 1;
    hipp = restrict(hipp, ExpKeys.TimeOnTrack, ExpKeys.TimeOffTrack);
end

%% remove dc shifts in voltage
if cfg_csc.detrend == 1;
    ofcdata = locdetrend(ofc.data, 1/ofc_dt, [1 0.5]);
    vstrdata = locdetrend(vstr.data, 1/vstr_dt, [1 0.5]);
    ofc.data = ofcdata'; % transpose after ofcdata because vdm lab tsd's are in the format [tvec = n x 1; data = 1 x n]
    vstr.data = vstrdata'; 
    if cfg_csc.hippflag == 1;
        hippdata = locdetrend(hipp.data, 1/hipp_dt, [1 0.5]);
        hipp.data = hippdata'; 
    end
end

%% First Order Differencing
if cfg_csc.diffdata == 1;   % removes much of the autocorrelation (but not all), if this is desired (for Granger analysis)
    ofc.tvec = ofc.tvec(1:end-1);
    ofc.data = diff(ofc.data);
    
    vstr.tvec = vstr.tvec(1:end-1);
    vstr.data = diff(vstr.data);
    
    if cfg_csc.hippflag == 1;
        hipp.tvec = hipp.tvec(1:end-1);
        hipp.data = diff(hipp.data);
    end
end

%% Checks
if sum(isnan(ofc.data))>0;
    warning('1 or more NaNs in OFC data');
end
if sum(isnan(vstr.data))>0;
    warning('1 or more NaNs in VSTR data');
end
if cfg_csc.hippflag == 1;
    if sum(isnan(hipp.data)) > 0;
        warning('1 or more NaNs in HIPP data');
    end
end

end

function [ofc, vstr, hipp, csc_csc] = prepCSCs_new(cfg_in)

%2014-09-10. JJS. Loads and prepares the ofc and vstr CSC files that are
%specified as 'good' in the keys file for analysis of granger causality for
%a single session.
%2020-04-06. JJS. Modifications to work w/ MvdMlab codeset.
cfg_def = [];
cfg_def.VoltageConvFactor = 10^6;
cfg_def.decimateFactor = 2; % Downsampling the data. Default here is 2 (from 2kHz to 1kHz). 
cfg_def.detrend = 1;   % for removing slow DC shifts in voltage
cfg_def.diffdata = 0;  % for removing autocorrelation in the time series
cfg_def.doRestrict = 1; % Restirct CSC to track time only. 
cfg_def.hippflag = 0;  % Do we want to Load a hippocampal CSC? (for sessions with electrode in fissure). 0 = no. 1 = yes.
% process_varargin(varargin);
[csc_csc] = ProcessConfig(cfg_def,cfg_in);

tic
LoadExpKeys;

%% Load the CSCs
cfg_ofc.fc = {ExpKeys.OFCcsc};
cfg_ofc.VoltageConvFactor = cfg_def.VoltageConvFactor;
cfg_ofc.decimateByFactor = csc_csc.decimateFactor;
[ofc] = LoadCSC(cfg_ofc);
ofc_dt = median(diff(ofc.tvec));

cfg_vstr.fc = {ExpKeys.VSTRcsc};
cfg_vstr.VoltageConvFactor = cfg_def.VoltageConvFactor;
cfg_vstr.decimateByFactor = csc_csc.decimateFactor;
[vstr] = LoadCSC(cfg_vstr);
vstr_dt = median(diff(vstr.tvec));

assert(length(ofc.data)==length(vstr.data));
assert(ofc_dt == vstr_dt);

if strcmp('NaN', ExpKeys.HIPPcsc) == 0 && csc_csc.hippflag == 1;
    cfg_hipp.fc = ExpKeys.HIPPcsc;
    cfg_hipp.VoltageConvFactor = cfg_def.VoltageConvFactor;
    cfg_hipp.decimateByFactor = csc_csc.decimateFactor;
    [hipp] = LoadCSC(cfg_hipp);
    hipp_dt = median(diff(hipp.tvec));
end

%% restrict to time on track
if cfg_def.doRestrict == 1;
    ofc = restrict(ofc, ExpKeys.TimeOnTrack, ExpKeys.TimeOffTrack);
    vstr = restrict(vstr, ExpKeys.TimeOnTrack, ExpKeys.TimeOffTrack);
end

if cfg_def.doRestrict == 1 && cfg_def.hippflag == 1;
    hipp = restrict(hipp, ExpKeys.TimeOnTrack, ExpKeys.TimeOffTrack);
end

%% remove dc shifts in voltage
if csc_csc.detrend == 1;
    ofcdata=locdetrend(ofc.data, 1/ofc_dt, [1 0.5]);
    vstrdata=locdetrend(vstr.data, 1/vstr_dt, [1 0.5]);
    ofc = tsd(ofc.tvec, ofcdata);
    vstr = tsd(vstr.tvec, vstrdata);
    if csc_csc.hippflag == 1;
        hippdata = locdetrend(hipp.data, 1/hipp_dt, [1 0.5]);
        hipp = tsd(hipp.tvec, hippdata);
    else
        hipp = [];
    end
end

%% First Order Differencing
if csc_csc.diffdata == 1;  
    ofc = tsd(ofc.tvec(1:end-1), diff(ofc.data));
    vstr = tsd(vstr.tvec(1:end-1), diff(vstr.data));
end

%% Checks
if sum(isnan(ofc.data))>0;
    warning('1 or more NaNs in OFC data');
end
if sum(isnan(vstr.data))>0;
    warning('1 or more NaNs in VSTR data');
end
if csc_csc.hippflag == 1;
    if sum(isnan(hipp.data))>0;
        warning('1 or more NaNs in HIPP data');
    end
end
disp('time to load CSC')
toc

end

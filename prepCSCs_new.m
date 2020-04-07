function [ofc, vstr, hipp, dt, hippflag, decimatefactor, detrend, diffdata] = prepCSCs_new
%2014-09-10. JJS. Loads and prepares the ofc and vstr CSC files that are
%specified as 'good' in the keys file for analysis of granger causality for
%a single session.
%2020-04-06. JJS. Modifications to work w/ MvdMlab codeset.

decimatefactor = 2; % decimates by taking every nth sample
filtorder = 8;
detrend = 1;
diffdata = 0;
doRestrict = 1;
% process_varargin(varargin);

hippflag = 0;
tic
EvalKeys;
% Load the CSCs
cfg_ofc.fc = {ExpKeys.OFCcsc};
cfg_ofc.VoltageConvFactor = 10^6;
cfg_ofc.decimateByFactor = 2;
[ofc] = LoadCSC(cfg_ofc);

cfg_vstr.fc = {ExpKeys.VSTRcsc}; 
cfg_vstr.VoltageConvFactor = 10^6;
cfg_vstr.decimateByFactor = 2;
[vstr] = LoadCSC(cfg_vstr); 

if strcmp('NaN', ExpKeys.HIPPcsc) == 0 && hippflag == 1;
    cfg_hipp.fc = ExpKeys.HIPPcsc;
    [hipp] = LoadCSC(cfg_hipp); 
end

% restrict to time on track
if doRestrict == 1;
    ofc = restrict(ofc, ExpKeys.TimeOnTrack, ExpKeys.TimeOffTrack);
    vstr = restrict(vstr, ExpKeys.TimeOnTrack, ExpKeys.TimeOffTrack);
end

if doRestrict == 1 && hippflag == 1;
    hipp = restrict(hipp, ExpKeys.TimeOnTrack, ExpKeys.TimeOffTrack);
end
assert(length(ofc.data)==length(vstr.data));
ofc_dt = median(diff(ofc.tvec)); 
vstr_dt = median(diff(vstr.tvec));
assert(ofc_dt == vstr_dt);

% remove dc shifts in voltage
if detrend == 1;
    ofcdata=locdetrend(ofc.data, 1/ofc.dt, [1 0.5]);
    vstrdata=locdetrend(vstr.data, 1/vstr.dt, [1 0.5]);
    ofc = tsd(ofc.range, ofcdata);
    vstr = tsd(vstr.range, vstrdata);
    if hippflag == 1;
        hippdata = locdetrend(hipp.data, 1/hipp.dt, [1 0.5]);
        hipp = tsd(hipp.range, hippdata);
    end
end

if diffdata == 1;
    % First Order Differencing
    ofc = tsd(ofc.T(1:end-1), diff(ofc.D));
    vstr = tsd(vstr.T(1:end-1), diff(vstr.D));
end

dt = 1/ofc.dt;
assert((1/ofc.dt)==(1/vstr.dt))

if sum(isnan(ofc.D))>0;
    warning('1 or more NaNs in OFC data');
end
if sum(isnan(vstr.D))>0;
    warning('1 or more NaNs in VSTR data');
end
if hippflag == 1;
    if sum(isnan(hipp.D))>0;
        warning('1 or more NaNs in HIPP data');
    end
end

toc

end

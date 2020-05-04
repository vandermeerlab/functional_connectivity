function sd = DDinit_new(varargin)

% sd = DDinit(fd)
%
% DelayDiscounting (DD) task initialization function
% checks and loads keys, video-tracker-1, video-tracker-2, mat file, events file, and spikes
%
% ADR 2011-12

%conditionals
restrictSessionToTrack = 0;
VT1 = 1;  %1 = load VT1, 0 = don't
VT2 = 1;  %1 = load VT2, 0 = don't
Spikes = 1;  %1 = load spikes, 0 = don't
limitSpikes = 1;
DD  = 1;  %1 = load *DD.mat, 0 = don't
Use__Ts = false;
if mod(length(varargin),2)~=0;
    varargin{end+1}=nan;
end
process_varargin(varargin);

fd = pwd;

assert(exist(fd, 'dir')==7, 'Cannot find directory %s.', fd);

[~, SSN, ~] = fileparts(fd);
sd.SSN = SSN;

% -----------------------
% KEYS
% -----------------------
keysfn = [strrep(SSN, '-', '_') '_keys'];
assert(exist(keysfn, 'file')==2, 'Cannot find keys file %s.', keysfn);
eval(keysfn);
sd.ExpKeys = ExpKeys;
sd.ExpKeys.SSN = SSN;
sd.ExpKeys.fd = fd;

assert(~iscell(ExpKeys.Behavior), 'Multiple Behaviors');

%------------------------
% VIDEO TRACKING
%------------------------
% Video-tracker-1
if VT1==1
    W = warning();
    warning off MATLAB:unknownObjectNowStruct
    vtfn = fullfile(fd, [SSN '-vt.mat']);
    assert(exist(vtfn, 'file')==2, 'Cannot find vt file %s.', vtfn);
    if exist(vtfn, 'file')
        load(vtfn);
        if exist('Vt', 'var'), x = Vt.x; y = Vt.y; end
%         if isstruct(x); x = tsd(x); end                   % TSD in vdml code changes x. 
%         if isstruct(y); y = tsd(y); end
%         [x,y] = LoadVT_lumrg(strcat(SSN,'-VT1.nvt'));
        if restrictSessionToTrack == 1;
            sd.x = restrict(x, ExpKeys.TimeOnTrack, ExpKeys.TimeOffTrack);
            sd.y = restrict(y, ExpKeys.TimeOnTrack, ExpKeys.TimeOffTrack);
        else
            sd.x = x;
            sd.y = y;
        end
    else
        warning('FileNotFound: No VT file found.');
    end
    warning(W);
end
% Video-tracker-2
if VT2==1
    W = warning();
    warning off MATLAB:unknownObjectNowStruct
    vtfn = fullfile(fd, [SSN '-vt2.mat']);
    assert(exist(vtfn, 'file')==2, 'Cannot find vt file %s.', vtfn);
    if exist(vtfn, 'file')
        load(vtfn);
        if exist('Vt', 'var'), x = Vt.x; y = Vt.y; end
%         if isstruct(x); x = tsd(x); end
%         if isstruct(y); y = tsd(y); end
        if restrictSessionToTrack == 1;
            sd.x2 = restrict(x, ExpKeys.TimeOnTrack, ExpKeys.TimeOffTrack);
            sd.y2 = restrict(y, ExpKeys.TimeOnTrack, ExpKeys.TimeOffTrack);
        else
            sd.x2 = x;
            sd.y2 = y;
        end
    else
        warning('FileNotFound: No VT2 file found.');
    end
    warning(W);
end

%-------------------------
% EVENTS
%-------------------------
eventsfn = fullfile(fd, [SSN '-events.Nev']);
if ~exist(eventsfn, 'file')==1; warning('No events file found'); end

% assert(exist(eventsfn, 'file')==2, 'Cannot find events file %s.', eventsfn);

%-------------------------
% SPIKES
%-------------------------
if Spikes ==1
    fc = FindFiles('*.t', 'CheckSubdirs', 0);
    if Use__Ts; fc = cat(1, fc, FindFiles('*._t', 'CheckSubdirs',0));  end %#ok<UNRCH>
    %     cfg_in.fc = {fc};
    S = LoadSpikes([]);
    %% THIS SECTION NEEDS TO BE REWORKED
%     L = zeros(1,length(S));
%     if restrictSessionToTrack == 1;
%         for iC = 1:length(S)
%             S.t{iC} = restrict(S.t{iC}, ExpKeys.TimeOnTrack, ExpKeys.TimeOffTrack);
%             L(iC) = length(S.t{iC}.data);
%         end
%     end
%     if limitSpikes == 1;
%         keep = L>100;
%         sd.S = S(keep);
%     else sd.S = S;
%     end
%     
%     sd.fc = fc(keep);
%     sd.fn = {};
%     for iC = 1:length(fc(keep))
%         [~,sd.fn{iC}] = fileparts(sd.fc{iC});
%         sd.fn{iC} = strrep(sd.fn{iC}, '_', '-');
%     end
%     sd.fn = sd.fn';
    sd.S = S;
end


%-------------------------
% TARGET ID
%-------------------------
% make array of length n cells with 1 to specify OFC and 2 to specify vStr
% Added by JJS on 3/28/2012
% NOTE: need to add a third option, '3', b/c we found out that for R214 I overshot ofc
% and ended up in piriform cortex. JJS. 2013-11-14.
if Spikes
    TargetID = nan(length(sd.S),1);
    for iS = 1:length(sd.S);
        [~, name, ~] = fileparts(fc{iS,1});
        TTposition = regexp(name, '-TT');
        if ~isempty(regexp(name, '-TT', 'once'));
            TTposition = regexp(name, '-TT');
        else TTposition = regexp(name, '-tt');
        end
        Tetrode = name(TTposition+3:TTposition+4);
        Tetrode = str2num(Tetrode);
        if sd.ExpKeys.TetrodeTargets(1,Tetrode)==1;
            TargetID(iS) = 1; %OFC
        elseif sd.ExpKeys.TetrodeTargets(1,Tetrode)==2;
            TargetID(iS) = 2; %vStr
        else
            SSNstr = GetSSN('SingleSession');
            fprintf('Unknown target.  Check keys file for session %s', SSNstr);
            error('Problem with sd.TargetID'); %#ok<LTARG>
        end
    end
    sd.OFC = sd.S(TargetID==1);
    sd.vStr = sd.S(TargetID==2);
    sd.TargetID = TargetID;
    sd.Target = cell(1,length(sd.S));
    for iT = 1:length(sd.Target);
        if sd.TargetID(iT) == 1;
            sd.Target{iT} = 'OFC';
        end
        if sd.TargetID(iT) == 2;
            sd.Target{iT} = 'Striatum';
        end
    end
end


%-----------------------
% *DD.mat File
%-----------------------

if DD==1
    DD = fullfile(fd, [SSN '-DD.mat']);
    assert(exist(DD, 'file')==2, 'Cannot find *DD.mat file %s.', DD);
    if exist(DD, 'file')
        load(DD);
        assert(length(ZoneIn)==TotalLaps); % Added by JJS. 2012-10-28
        %process *DD.mat file
        %         [sd.FeedersFired, FeederTimes, skip, ZoneDelay] = DDFixFeedersFired(FeederTimes, FeedersFired, ZoneIn, TotalLaps, ZoneDelay);
%         [FeedersFired, FeederTimes] = DD6_FixFeedersFired(FeederTimes, FeedersFired, ZoneIn, TotalLaps); %#ok<NODEF>
        
        % Edit. 2012-08-06 JJS. Fixes problem with zeros in FeedersFired (instead of NaNs) for laps in which the rat skips the feeder fire
        lapskips = FeedersFired==0;
        FeedersFired(lapskips==1)= NaN;
        sd.FeedersFired = FeedersFired;
        ZoneDelay(lapskips==1) = NaN;
        FeederTimes(lapskips==1)=NaN;
        
        [sd.DelayZone, ~,~] = DD_getWorld(World);
        if sd.DelayZone == 4; nonDelayZone = 3; end
        if sd.DelayZone == 3; nonDelayZone = 4; end
        sd.nonDelayZone = nonDelayZone;
        %convert [us] to [s]
        sd.FeederTimes = FeederTimes*(10^-6);
        sd.EnteringZoneTime = EnteringZoneTime*(10^-6);
        sd.ExitZoneTime = ExitZoneTime*(10^-6);
        
        sd.ZoneDelay = ZoneDelay;
        % Added 2012-05-30 by JJS
        AdjDel = sd.ZoneDelay;
        for iL = 1:length(AdjDel)
            if AdjDel(1)==1;
                AdjDel(1) = NaN;
            end
            if AdjDel(iL)==1;
                AdjDel(iL) = AdjDel(iL-1);
            else
                AdjDel(iL) = AdjDel(iL);
            end
        end
        sd.ZoneIn = ZoneIn;
        
        sd.AdjDel_all_laps = AdjDel;
        sd.AdjDel_del_laps = sd.ZoneDelay(sd.ZoneIn==sd.DelayZone);
        %         sd.AdjDel_del_laps = sd.AdjDel_del_laps(~isnan(sd.AdjDel_del_laps));
        
        sd.World = World;
        sd.TotalLaps = TotalLaps;
        
        % Added by JJS 2013-09-16. Gives the ave. final adjust. delay.
        lastTwenty = sd.TotalLaps - 20;
        sd.FinalAveDelay = nanmean(sd.AdjDel_all_laps(lastTwenty:end));
        %         sd.FL = sd.FeederTimes(sd.FeedersFired == 3);
        %         sd.FR = sd.FeederTimes(sd.FeedersFired == 1);
        % As defined here, sd.FL & sd.FR (as well as sd.D & sd.nD) exclude
        % laps without feeder fire events(i.e. skips or backtracking or
        % video errors). They do not contain NaNs. sd.FL + sd.FR may not
        % equal sd.TotalLaps.
        
        % Changed by JJS on 2012-10-28. Commented out the above formula in
        % place of the one below. Here, sd.FL and sd.FR can contain NaNs,
        % and sd.FL + sd.FR = sd.TotalLaps in all cases.
        % Note: confirm that sd.ZoneIn is always equal to sd.TotalLaps.
        sd.FL = sd.FeederTimes(sd.ZoneIn == 3);
        sd.FR = sd.FeederTimes(sd.ZoneIn == 4);
        
        
        
        if exist('EnteringCPTime','var'); sd.EnteringCPTime=EnteringCPTime; end;
        if exist('ExitingCPTime','var'); sd.ExitingCPTime=ExitingCPTime; end;
        if exist('FeederSkip','var'); sd.FeederSkip=FeederSkip; end;
        if exist('BackwardsTimes','var'); sd.BackwardsTimes=BackwardsTimes; end;
        if exist('EnteringSoMTime','var'); sd.EnteringSoMTime=EnteringSoMTime;end;
        if exist('ExitingSoMTime','var'); sd.ExitingSoMTime=ExitingSoMTime; end;
        
        
    else
        warning('FileNotFound: No *DD.mat file found.');
    end
    
    if sd.DelayZone == 3; sd.DelaySide = 'Left'; end
    if sd.DelayZone == 4; sd.DelaySide = 'Right'; end
    
end
% Added by JJS 2012-05-02

sd.EnteringZoneTimeFL = sd.EnteringZoneTime(sd.ZoneIn==3);
sd.EnteringZoneTimeFR = sd.EnteringZoneTime(sd.ZoneIn==4);

if sd.World.nPleft == 3;
    sd.delayedSide = 'Left';
    sd.nondelayedSide = 'Right';
else
    sd.delayedSide = 'Right';
    sd.nondelayedSide = 'Left';
end

% Added by JJS 2012-05-31
for iL = 1:sd.TotalLaps;
    if sd.ZoneIn(iL) == sd.DelayZone;
        Ltype(iL) = 1; % 1 equals a DELAY side lap
    end
    if sd.ZoneIn(iL) == sd.nonDelayZone;
        Ltype(iL) = 0; % 0 equals non-DELAY side lap
    end
end
sd.FeederType = Ltype;

% Added by JJS 2012-06-01
% Feeder times for the delayed and non-delayed sides.
if sd.DelayZone == 4;
    sd.D = sd.FR; sd.nD = sd.FL; end
if sd.DelayZone == 3;
    sd.D = sd.FL; sd.nD = sd.FR; end

% Added by JJS 2012-07-18
if Spikes
    sd.RunTime = sd.x.t(end) - sd.x.t(1);
    %     for iC = 1:length(sd.S);
    %         sd.aveFR{iC} = length(sd.S{iC}.T)/sd.RunTime;
    %     end
    %     if isempty(sd.S);
    %         sd.aveFR = NaN;
    %     end
    %
    %     for iC = 1:length(sd.OFC);
    %         sd.aveFR_OFC{iC} = length(sd.OFC{iC}.T)/sd.RunTime;
    %     end
    %     if isempty(sd.OFC);
    %         sd.aveFR_OFC = NaN;
    %     end
    %
    %     for iC = 1:length(sd.vStr);
    %         sd.aveFR_vStr{iC} = length(sd.vStr{iC}.T)/sd.RunTime;
    %     end
    %     if isempty(sd.vStr);
    %         sd.aveFR_vStr = NaN;
    %     end
    % Added by JJS 2012-08-02
%     SpikeCount = barspikes(sd.S);    % commented out 4/2020
%     sd.aveFR= SpikeCount / sd.RunTime;
end

% Added by JJS 2012-08-22
% Sort the tracking data in case timestamps are out of order.
% sd.x = tsd(sort(sd.x.t),sd.x.t);  % commented out 4/2020
% sd.y = tsd(sort(sd.y.t),sd.y.t);

% 208-02-15. Just trying to get this to work to look at some old stuff.
% sd.x2 = tsd(sort(sd.x2.range),sd.x2.data);
% sd.y2 = tsd(sort(sd.y2.range),sd.y2.data);

assert(length(sd.FeedersFired) == length(sd.FeederTimes));
% Added by JJS 2014-01-14
%  2014-02-02. JJS. Commented out. This categorize function needs work.
% if ~isempty(sd.vStr);
%     [MSP, TAN, PVB] = CategorizeStriatumFromDave(sd.vStr);
%     sd.MSP = MSP;
%     sd.TAN = TAN;
%     sd.FSI = PVB; % I prefer to call the putative PVB cells FSIs.
% else
%     sd.MSP = nan;
%     sd.TAN = nan;
%     sd.FSI = nan;
% end

% Added by JJS. 2013-05-08
% stored values for the estimated discount rates for each of my 6 rats.
% These estimates use the average final delay over all sessions and k =
% 2/d.
if strcmp(SSN(1:4),'R206'); sd.k = 0.3609; end
if strcmp(SSN(1:4),'R214'); sd.k = 0.4286; end
if strcmp(SSN(1:4),'R224'); sd.k = 0.3515; end
if strcmp(SSN(1:4),'R226'); sd.k = 0.3314; end
if strcmp(SSN(1:4),'R235'); sd.k = 0.2844; end
if strcmp(SSN(1:4),'R244'); sd.k = 0.3791; end

% Added by JJS. 2013-11-11.
for iT = 1:12;
    if sd.ExpKeys.TetrodeTargets(iT) == 1;
        sd.CSCtype{iT} = 'OFC';
    end
    if sd.ExpKeys.TetrodeTargets(iT) == 2;
        sd.CSCtype{iT} = 'Striatum';
    end
end

% Added by JJS. 2015-03-26.
if sd.DelayZone == 3; sd.Left = 'Delayed Side'; sd.Right = 'non-Delayed Side';
elseif sd.DelayZone == 4; sd.Left = 'non-Delayed Side'; sd.Right = 'Delayed Side';
else error('problem with zone ID')
end

% Added by JJS. 2015-04-16.

load(strcat(SSN, '-CurvatureTimes.mat'));
VTEindex = ~isnan(eventStarts(:,1));
nonVTEindex = isnan(eventStarts(:,1));
sd.VTElaps = find(VTEindex==1);
sd.nonVTElaps = find(VTEindex==0);
assert(length(sd.VTElaps)+length(sd.nonVTElaps)==sd.TotalLaps);

sd.VTEtimes = eventStarts(VTEindex,1);
sd.nonVTEtimes = eventStarts(nonVTEindex,1);
assert(length(sd.VTEtimes)+length(sd.nonVTEtimes)==sd.TotalLaps);

% function [ sd ] = DDGetCPpasses( sd, varargin )
% 2012-04-26 AndyP
% get choice point entering and exiting times
% Added by JJS. 2015-04-16.
splitTime=5;
X=sd.x.data;
Y=sd.y.data;
time=sd.x.t;
L0 = [sd.x.t(1), sd.ExitZoneTime];
L1 = [sd.EnteringZoneTime, sd.x.t(end)];
keep = 195 < X & X < 205 & 160 < Y & Y < 240;
if ~isempty(keep); [Lstart,~] = FindAnyLap([sd.x.t(1); time(keep)], X(keep),'splitTime', splitTime); end
EnteringCPTime = nan(1,sd.TotalLaps);
ExitingCPTime = nan(1,sd.TotalLaps);
if ~isempty(keep);
    for iL = 1:sd.TotalLaps;
        LstartToUse = find(Lstart > L0(iL) & Lstart < L1(iL));
        if isempty(LstartToUse)
        else
            if length(LstartToUse)>1; LstartToUse=LstartToUse(end); end
            EnteringCPTime(1,iL) =Lstart(LstartToUse);
            ExitingCPTime(1,iL) = L1(iL);
            assert(ExitingCPTime(iL)-EnteringCPTime(iL)>0,'negative times');
        end
    end
else error('no tracking data');
end
sd.EnteringCPTime = EnteringCPTime;
sd.ExitingCPTime = ExitingCPTime;
sd.CPpassDuration = sd.ExitingCPTime - sd. EnteringCPTime;

%% Added by JJS. 2013-03-24. Code is from AEP's DDinit function. It determines the coordinates for the feeder zones.
% * double check to make sure Andy's coordinates line up with
% (approximately) with the coordinates from my recording sessions.
%-----------------------
% DD.mat File
%-----------------------
% if DD==1
DD = fullfile(fd, [SSN '-DD.mat']);
assert(exist(DD, 'file')==2, 'Cannot find *DD.mat file %s.', DD);
if exist(DD, 'file')
    load(DD); % process *DD.mat file
    %
    % 			[FeedersFired, FeederTimes] = DD6_FixFeedersFired(Fee		if ~strcmp(sd.ExpKeys.Behavior,'FPT Delay Discounting');derTimes, FeedersFired, ZoneIn, TotalLaps); %#ok<NODEF>
    if datenum(SSN(6:end),'yyyy-mm-dd')<datenum('2011-09-05','yyyy-mm-dd')
        %fix bug in feeders fired and feeder times 2012-03-17 AndyP
        sd.Coord.SoM_x = 280; % Start of Maze <x,y>
        sd.Coord.SoM_y = 209;
        sd.Coord.CP_x = 141;  % Choice point <x,y>
        sd.Coord.CP_y = 209;
        sd.Coord.LF_x = 141;  % Left feeder <x,y>
        sd.Coord.LF_y = 337;
        sd.Coord.RF_x = 141;  % Right Feeder <x,y>
        sd.Coord.RF_y = 81;
        sd.InZoneDistance = [77 77 77 77]; %radius of zone [1 2 3 4] (# pixels)
    elseif datenum(SSN(6:end),'yyyy-mm-dd')>datenum('2011-09-05','yyyy-mm-dd')
        sd.Coord.SoM_x = 280; % Start of Maze <x,y>
        sd.Coord.SoM_y = 209;
        sd.Coord.CP_x = 141;  % Choice point <x,y>
        sd.Coord.CP_y = 209;
        sd.Coord.LF_x = 141;  % Left feeder <x,y>
        sd.Coord.LF_y = 337;
        sd.Coord.RF_x = 141;  % Right Feeder <x,y>
        sd.Coord.RF_y = 81;
        sd.InZoneDistance = [60 60 60 60]; %radius of zone [1 2 3 4] (# pixels)
    end
    % 		else
    % 		end
    %--------------
    % 		[sd.DelayZone, ~,~] = DD_getWorld(World);
    % 		sd.ZoneDelay = ZoneDelay;
    % 		sd.ZoneIn = ZoneIn;
    % 		sd.World = World;
    % 		sd.TotalLaps = TotalLaps;
    %--------------
    % 		if exist('EnteringZoneTime','var'); sd.EnteringZoneTime = EnteringZoneTime*(10^-6); end
    % 		if exist('ExitZoneTime','var'); sd.ExitZoneTime = ExitZoneTime*(10^-6); end
    % 		if exist('FeedersFired','var'); sd.FeedersFired=FeedersFired; end
    % 		if exist('FeederTimes','var'); sd.FeederTimes = FeederTimes*(10^-6); end
    % 		if exist('BackwardsTimes','var'); sd.BackwardsTimes=BackwardsTimes*(10^-6); end
    % 		if exist('EnteringSoMTime','var'); sd.EnteringSoMTime=EnteringSoMTime*(10^-6); end
    % 		if exist('ExitingSoMTime','var'); sd.ExitingSoMTime=ExitingSoMTime*(10^-6); end
    % 		if exist('EnteringCPTime','var'); sd.EnteringCPTime=EnteringCPTime.*(10^-6); end
    % 		if exist('ExitingCPTime','var'); sd.ExitingCPTime=ExitingCPTime.*(10^-6); end
    %
    % 		if exist('FeederSkip','var'); sd.FeederSkip=FeederSkip; end
    % 		if exist('SwitchLap','var'); sd.SwitchLap=SwitchLap; end
    
else
    warning('FileNotFound: No *DD.mat file found.');
end
% end
%--------------


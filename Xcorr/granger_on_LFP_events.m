function [f12, f21, params, fvec, pval, sig]  = granger_on_LFP_events(X, varargin)
% 2018-09-08. JJS. 'granger_ver2'
% Revised 2018.

% Issues: AIC and BIC decrease monotonically. Best model order is always the one with the max number of lags.
%
%% granger causality
% INPUTS

% OUTPUTS
% 1 is OFC and 2 is vSTR

% F12 - time series GC for ofc -> vstr
% F21 - time series GC for vstr -> ofc
% f12 (1 x nfreqs) - spectral GC for ofc -> vstr
% f21 (1 x nfreqs) - spectral GC for vstr -> ofc

% fvec - the array of frequencies used for spectral GC calculation
% fs - sampling rate, in Hertz
% params - structure with the parameters used to run the GC calculation
% pval - pvalue for significance of time series GC (versus the null hypothesis)
% sig - 0 or 1. 0 = not sig at p = .05. 1 = sig at p = .05.

% VARARGIN OPTIONS
% see parameters below
%% Functions used in this function
% tsdata_to_infocrit(D,momax,icregmode);
% tsdata_to_var(D,morder,regmode);
% var_to_autocov(A,SIG,acmaxlags);
% autocov_to_spwcgc(G,fres);
%% Parameters
doPlotSig = 0;
doPlot = 1;

LoadExpKeys;
% ofc = ofc.restrict(ExpKeys.TimeOnTrack, ExpKeys.TimeOffTrack);
% vstr = vstr.restrict(ExpKeys.TimeOnTrack, ExpKeys.TimeOffTrack);
% [ofc, vstr, ~, ~, ~, ~, ~, ~] = prepCSCs;
% D = ofc.data; D = D';
% D(2,:) = vstr.data;

D = X.eventData; 
params = [];
nvars = size(D,1);
ntrials = size(D, 3);     % number of gamma events
nobs    = size(D,2);   % number of observations per trial
% regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
regmode   = [];  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)
morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 20;     % maximum model order for model order estimation
acmaxlags = [];   % maximum autocovariance lags (empty for automatic calculation)
tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')
% fs        = 1/ofc.dt;    % sample rate (Hz)
% assert(1/ofc.dt == 1/vstr.dt);
fres = 300;
fs = X.eventStats.usr.fs;

process_varargin(varargin);

%% Parameters
params.regmode = regmode;
params.icregmode = icregmode;
params.morder = morder;
params.momax = momax;
params.acmaxlags = acmaxlags;
params.tstat = tstat;
params.alpha = alpha;
params.mhtc = mhtc;
params.ntrials = ntrials;
params.nvars = nvars;
params.fres = fres;
params.fs = fs; % sampling rate

%% Model order estimation (A2)
% Calculate information criteria up to specified maximum model order.
ptic('\n*** tsdata_to_infocrit\n');
[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(D,momax,icregmode);
% disp('AIC'); disp(AIC);figure
% disp('BIC'); disp(BIC);
ptoc('*** tsdata_to_infocrit took ');
% Plot information criteria.
if doPlot ==1;
    figure(1); 
    plot_tsdata([AIC BIC]',{'AIC','BIC'},1/fs);
    title('Model order estimation');
end

amo = size(D,1); % actual model order

fprintf('\nbest model order (AIC) = %d\n',moAIC);
fprintf('best model order (BIC) = %d\n',moBIC);
fprintf('actual model order     = %d\n',amo);
% Select model order.
if     strcmpi(morder,'actual')
    morder = amo;
    fprintf('\nusing actual model order = %d\n',morder);
elseif strcmpi(morder,'AIC')
    morder = moAIC;
    fprintf('\nusing AIC best model order = %d\n',morder);
elseif strcmpi(morder,'BIC')
    morder = moBIC;
    fprintf('\nusing BIC best model order = %d\n',morder);
else
    fprintf('\nusing specified model order = %d\n',morder);
end

%% VAR model estimation (A2)
% Estimate VAR model of selected order from data.
ptic('\n*** tsdata_to_var... ');
[A,SIG] = tsdata_to_var(D,morder,regmode);
ptoc;
% Check for failed regression
assert(~isbad(A),'VAR estimation failed');
% NOTE: at this point we have a model and are finished with the data! - all
% subsequent calculations work from the estimated VAR parameters A and SIG.

%% Autocovariance calculation (A5)

% The autocovariance sequence drives many Granger causality calculations (see
% next section). Now we calculate the autocovariance sequence G according to the
% VAR model, to as many lags as it takes to decay to below the numerical
% tolerance level, or to acmaxlags lags if specified (i.e. non-empty).
ptic('*** var_to_autocov... ');
[G,info] = var_to_autocov(A,SIG,acmaxlags);
ptoc;
% The above routine does a LOT of error checking and issues useful diagnostics.
% If there are problems with your data (e.g. non-stationarity, colinearity,
% etc.) there's a good chance it'll show up at this point - and the diagnostics
% may supply useful information as to what went wrong. It is thus essential to
% report and check for errors here.

var_info(info,true); % report results (and bail out on error)

%% Granger causality calculation: frequency domain (A14)

% Calculate spectral pairwise-conditional causalities at given frequency
% resolution - again, this only requires the autocovariance sequence.
ptic('\n*** autocov_to_spwcgc... ');
f = autocov_to_spwcgc(G,fres);
ptoc;
% Check for failed spectral GC calculation
assert(~isbad(f,false),'spectral GC calculation failed');
figure(2)
fvec = plot_spw(f,fs); close;

f21 = squeeze(f(1,2,:));
f12 = squeeze(f(2,1,:)); 

% Plot spectral causal graph.
if doPlot == 1;
    figure(3); clf;
    plot_spw(f,fs);
    subplot(2,2,2)
    c = axis;
    axis([c(1) 100 c(3) c(4)]);
    title('Spectral pairwise G-causality: vstr --> ofc')
    subplot(2,2,3)
    c = axis;
    axis([c(1) 100 c(3) c(4)]);
    title('Spectral pairwise G-causality: ofc --> vstr')
end
%% Granger causality calculation: time domain  (<mvgc_schema.html#3 |A13|>)

% Calculate time-domain pairwise-conditional causalities - this just requires
% the autocovariance sequence.

ptic('*** autocov_to_pwcgc... ');
F = autocov_to_pwcgc(G);
ptoc;

disp(F);
% Check for failed GC calculation

assert(~isbad(F,false),'GC calculation failed');

% Significance test using theoretical null distribution, adjusting for multiple hypotheses.

pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2,tstat); % take careful note of arguments!
sig  = significance(pval,alpha,mhtc);

disp(pval)
disp(sig)
% Plot time-domain causal graph, p-values and significance.

if doPlotSig == 1;
    figure(3); clf;
    subplot(1,3,1);
    plot_pw(F);
    title('Pairwise-conditional GC');
    subplot(1,3,2);
    plot_pw(pval);
    title('p-values');
    subplot(1,3,3);
    plot_pw(sig);
    title(['Significant at p = ' num2str(alpha)])
end

F12 = F(2,1); disp(strcat('ofc-->vstr =', num2str(F12)));
F21 = F(1,2); disp(strcat('vstr-->ofc =', num2str(F21)));

% For good measure we calculate Seth's causal density (cd) measure - the mean
% pairwise-conditional causality. We don't have a theoretical sampling
% distribution for this.

cd = mean(F(~isnan(F)));

fprintf('\ncausal density = %f\n',cd);

%% Granger causality calculation: frequency domain -> time-domain  (<mvgc_schema.html#3 |A15|>)

% Check that spectral causalities average (integrate) to time-domain
% causalities, as they should according to theory.

fprintf('\nchecking that frequency-domain GC integrates to time-domain GC... \n');
Fint = smvgc_to_mvgc(f); % integrate spectral MVGCs
mad = maxabs(F-Fint);
madthreshold = 1e-5;
if mad < madthreshold
    fprintf('maximum absolute difference OK: = %.2e (< %.2e)\n',mad,madthreshold);
else
    fprintf(2,'WARNING: high maximum absolute difference = %e.2 (> %.2e)\n',mad,madthreshold);
end

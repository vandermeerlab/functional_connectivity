function [CROSSCORRordered] = sort_crosscorr_single_session(X, keeplist, varargin)

%   inputs:     X.crosscorr - an n session by m X.lags matrix of X.crosscorr values

%   outputs:    cc_data_sorted - same as above, but sorted top to bottom
%   with max X.lags increasing from left to right

doPlot = 1;
FontSize = 15;
process_varargin(varargin);

if ~isempty(keeplist);
    X.max_crosscorr_lag = X.max_crosscorr_lag(X.keep & keeplist)';
    X.crosscorr = X.crosscorr(:,X.keep & keeplist)';
% elseif isstring(keeplist) & strcmp(keeplist, 'All') == 1; 
else
    X.max_crosscorr_lag = X.max_crosscorr_lag(X.keep)';
    X.crosscorr = X.crosscorr(:,X.keep)';
end

noNaNindex = ~isnan(X.max_crosscorr_lag);
temp = find(noNaNindex, 1, 'first');
X.lags = X.lags(temp,:);
CROSSCORRto_use = X.crosscorr(noNaNindex,:);
normCROSSCORR = nan(size(CROSSCORRto_use));
for iEvent = 1:size(CROSSCORRto_use,1);
    data = CROSSCORRto_use(iEvent,:);
    result = -1 + 2.*(data - min(data))./(max(data) - min(data)); % *from MATLAB central.
    normCROSSCORR(iEvent,:) = result;
    [C,I(iEvent)] = max(result);
    assert(C==1)
end
[B,IX] = sort(I);
for iEvent = 1:size(CROSSCORRto_use,1);
    maxLag(iEvent) = X.lags(B(iEvent));
end

CROSSCORRordered = nan(size(normCROSSCORR));
for iEvent = 1:size(CROSSCORRto_use,1);
    CROSSCORRordered(iEvent,:) = normCROSSCORR(IX(iEvent),:);
end
if doPlot == 1;
    imagesc(X.lags(1,:),1:size(CROSSCORRto_use,1),CROSSCORRordered)
    set(gca, 'FontSize', FontSize)
    xlabel('Lag (ms)', 'FontSize', FontSize)
    ylabel('Event Count', 'FontSize', FontSize)
%     title(sprintf('%s', num2str(X.cfg.low_freq),'-',num2str(X.cfg.high_freq), ' (Hz)'))
    c = axis; hold on
    line([0 0],[c(3) c(4)], 'color', 'k','LineWidth',3)
    for iEvent = 1:size(CROSSCORRto_use,1);
        plot(maxLag(iEvent),iEvent, '*', 'color', 'white')
    end
end

CROSSCORRordered = CROSSCORRordered';

% * http://www.mathworks.com/matlabcentral/answers/154075-how-to-scale-normalize-values-in-a-matrix-to-be-between-1-and-1

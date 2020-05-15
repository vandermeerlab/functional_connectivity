function Xcorr_subplot(X, varargin)

nBins = 10;
FontSize = 15;
process_varargin(varargin); 
clf; hold on
% SSN = GetSSN('SingleSession');
% sgtitle(strcat(SSN, ' : ', num2str(X.cfg.low_freq), '-', num2str(X.cfg.high_freq), 'Hz'));

colormap(jet)
% load(strcat(SSN,'-XCorr-lowgamma.mat'));
%% Crosscorr sorted heatmap of events
subplot(3,7,1)
sort_crosscorr_single_session(X, X.PreRun_LFPevents, 'doPlot', 1);
title('Pre-Run', 'FontSize', FontSize)
subplot(3,7,2)
sort_crosscorr_single_session(X, X.VTE_LFPevents, 'doPlot', 1); % VTE passes 
title('VTE passes', 'FontSize', FontSize)
subplot(3,7,3)
sort_crosscorr_single_session(X, X.nonVTE_LFPevents, 'doPlot', 1); % nonVTE passes 
title('nonVTE passes', 'FontSize', FontSize)
subplot(3,7,4)
sort_crosscorr_single_session(X, X.Approach_LFPevents, 'doPlot', 1);
title('Feeder Approach', 'FontSize', FontSize)
subplot(3,7,5)
sort_crosscorr_single_session(X, X.Waiting_LFPevents, 'doPlot', 1);
title('Reward Anticipation', 'FontSize', FontSize)
subplot(3,7,6)
sort_crosscorr_single_session(X, X.Reward_LFPevents, 'doPlot', 1);
title('Reward Consumption', 'FontSize', FontSize)
subplot(3,7,7)
sort_crosscorr_single_session(X, X.PostRun_LFPevents, 'doPlot', 1);
title('Post-run', 'FontSize', FontSize)

%% Average Crosscorr as a function of lag 
subplot(3,7,8)
shadedErrorBar(X.lags', nanmean(X.crosscorr(:,X.PreRun_LFPevents), 2), nanstderr(X.crosscorr(:,X.PreRun_LFPevents)'))
axis([-100 100 -.05 1.05]);
set(gca, 'FontSize', FontSize)
xlabel('lags (ms)', 'FontSize', FontSize)
ylabel('Xcorr','FontSize',FontSize)
[~, i] = max(nanmean(X.crosscorr(:,X.PreRun_LFPevents),2)); maxlag = X.lags(i); 
line([0 0], [-.05 1.05], 'color', 'k', 'LineWidth', 1, 'LineStyle', '-');
line([maxlag maxlag], [-.05 1.05], 'color', 'k', 'LineWidth', 1, 'LineStyle', '--');
text(-95,.95, strcat('maxlag=', num2str(maxlag, '%.1f')), 'color', 'k', 'FontSize', 8)

subplot(3,7,9)
shadedErrorBar(X.lags', nanmean(X.crosscorr(:,X.VTE_LFPevents), 2), nanstderr(X.crosscorr(:,X.VTE_LFPevents)'))
axis([-100 100 -.05 1.05]);
set(gca, 'FontSize', FontSize)
xlabel('lags (ms)', 'FontSize', FontSize)
ylabel('Xcorr','FontSize',FontSize)
[~, i] = max(nanmean(X.crosscorr(:,X.VTE_LFPevents),2)); maxlag = X.lags(i);
line([0 0], [-.05 1.05], 'color', 'k', 'LineWidth', 1, 'LineStyle', '-');
line([maxlag maxlag], [-.05 1.05], 'color', 'k', 'LineWidth', 1, 'LineStyle', '--');
text(-95,.95, strcat('maxlag=', num2str(maxlag, '%.1f')), 'color', 'k', 'FontSize', 8)

subplot(3,7,10)
shadedErrorBar(X.lags', nanmean(X.crosscorr(:,X.nonVTE_LFPevents), 2), nanstderr(X.crosscorr(:,X.nonVTE_LFPevents)'))
axis([-100 100 -.05 1.05]);
set(gca, 'FontSize', FontSize)
xlabel('lags (ms)', 'FontSize', FontSize)
ylabel('Xcorr','FontSize',FontSize)
[~, i] = max(nanmean(X.crosscorr(:,X.nonVTE_LFPevents),2)); maxlag = X.lags(i);
line([0 0], [-.05 1.05], 'color', 'k', 'LineWidth', 1, 'LineStyle', '-');
line([maxlag maxlag], [-.05 1.05], 'color', 'k', 'LineWidth', 1, 'LineStyle', '--');
text(-95,.95, strcat('maxlag=', num2str(maxlag, '%.1f')), 'color', 'k', 'FontSize', 8)

subplot(3,7,11)
shadedErrorBar(X.lags', nanmean(X.crosscorr(:,X.Approach_LFPevents), 2), nanstderr(X.crosscorr(:,X.Approach_LFPevents)'))
axis([-100 100 -.05 1.05]);
set(gca, 'FontSize', FontSize)
xlabel('lags (ms)', 'FontSize', FontSize)
ylabel('Xcorr','FontSize',FontSize)
[~, i] = max(nanmean(X.crosscorr(:,X.Approach_LFPevents),2)); maxlag = X.lags(i);
line([0 0], [-.05 1.05], 'color', 'k', 'LineWidth', 1, 'LineStyle', '-');
line([maxlag maxlag], [-.05 1.05], 'color', 'k', 'LineWidth', 1, 'LineStyle', '--');
text(-95,.95, strcat('maxlag=', num2str(maxlag, '%.1f')), 'color', 'k', 'FontSize', 8)

subplot(3,7,12)
shadedErrorBar(X.lags', nanmean(X.crosscorr(:,X.Waiting_LFPevents), 2), nanstderr(X.crosscorr(:,X.Waiting_LFPevents)'))
axis([-100 100 -.05 1.05]);
set(gca, 'FontSize', FontSize)
xlabel('lags (ms)', 'FontSize', FontSize)
ylabel('Xcorr','FontSize',FontSize)
[~, i] = max(nanmean(X.crosscorr(:,X.Waiting_LFPevents),2)); maxlag = X.lags(i);
line([0 0], [-.05 1.05], 'color', 'k', 'LineWidth', 1, 'LineStyle', '-');
line([maxlag maxlag], [-.05 1.05], 'color', 'k', 'LineWidth', 1, 'LineStyle', '--');
text(-95,.95, strcat('maxlag=', num2str(maxlag, '%.1f')), 'color', 'k', 'FontSize', 8)

subplot(3,7,13)
shadedErrorBar(X.lags', nanmean(X.crosscorr(:,X.Reward_LFPevents), 2), nanstderr(X.crosscorr(:,X.Reward_LFPevents)'))
axis([-100 100 -.05 1.05]);
set(gca, 'FontSize', FontSize)
xlabel('lags (ms)', 'FontSize', FontSize)
ylabel('Xcorr','FontSize',FontSize)
[~, i] = max(nanmean(X.crosscorr(:,X.Reward_LFPevents),2)); maxlag = X.lags(i);
line([0 0], [-.05 1.05], 'color', 'k', 'LineWidth', 1, 'LineStyle', '-');
line([maxlag maxlag], [-.05 1.05], 'color', 'k', 'LineWidth', 1, 'LineStyle', '--');
text(-95,.95, strcat('maxlag=', num2str(maxlag, '%.1f')), 'color', 'k', 'FontSize', 8)

subplot(3,7,14)
shadedErrorBar(X.lags', nanmean(X.crosscorr(:,X.PostRun_LFPevents), 2), nanstderr(X.crosscorr(:,X.PostRun_LFPevents)'))
axis([-100 100 -.05 1.05]);
set(gca, 'FontSize', FontSize)
xlabel('lags (ms)', 'FontSize', FontSize)
ylabel('Xcorr','FontSize',FontSize)
[~, i] = max(nanmean(X.crosscorr(:,X.PostRun_LFPevents),2)); maxlag = X.lags(i);
line([0 0], [-.05 1.05], 'color', 'k', 'LineWidth', 1, 'LineStyle', '-');
line([maxlag maxlag], [-.05 1.05], 'color', 'k', 'LineWidth', 1, 'LineStyle', '--');
text(-95,.95, strcat('maxlag= ', num2str(maxlag, '%.1f')), 'color', 'k', 'FontSize', 8)

%% Histogram of Max Lags 
subplot(3,7,15)
hist(X.max_crosscorr_lag(:,X.PreRun_LFPevents), nBins)
c = axis;
hold on
mean_h = mean(X.max_crosscorr_lag(:,X.PreRun_LFPevents)); 
median_h = median(X.max_crosscorr_lag(:,X.PreRun_LFPevents));
plot(mean_h, c(4), 'kV')
plot(median_h, c(4), 'rV')
xlabel('lags (ms)', 'FontSize', FontSize)
ylabel('Num. of Events', 'FontSize', FontSize)
set(gca, 'FontSize', FontSize)
axis([-100 100 0 c(4)]);
line([0 0], [0 c(4)], 'color', 'k', 'LineWidth', 2, 'LineStyle', '-');
text(.025,.95, strcat('mean= ', num2str(mean_h, '%.1f')), 'color', 'k', 'FontSize', 8, 'Units', 'Normalized')
text(.025,.90, strcat('med.= ', num2str(median_h, '%.1f')), 'color', 'k', 'FontSize', 8, 'Units', 'Normalized')

subplot(3,7,16)
hist(X.max_crosscorr_lag(:,X.VTE_LFPevents), nBins)
c = axis;
hold on
mean_h = mean(X.max_crosscorr_lag(:,X.VTE_LFPevents)); 
median_h = median(X.max_crosscorr_lag(:,X.VTE_LFPevents));
plot(mean_h, c(4), 'kV')
plot(median_h, c(4), 'rV')
xlabel('lags (ms)', 'FontSize', FontSize)
ylabel('Num. of Events', 'FontSize', FontSize)
set(gca, 'FontSize', FontSize)
axis([-100 100 0 c(4)]);
line([0 0], [0 c(4)], 'color', 'k', 'LineWidth', 2, 'LineStyle', '-');
text(.025,.95, strcat('mean= ', num2str(mean_h, '%.1f')), 'color', 'k', 'FontSize', 8, 'Units', 'Normalized')
text(.025,.90, strcat('med.= ', num2str(median_h, '%.1f')), 'color', 'k', 'FontSize', 8, 'Units', 'Normalized')

subplot(3,7,17)
hist(X.max_crosscorr_lag(:,X.nonVTE_LFPevents), nBins)
c = axis;
hold on
mean_h = mean(X.max_crosscorr_lag(:,X.nonVTE_LFPevents)); 
median_h = median(X.max_crosscorr_lag(:,X.nonVTE_LFPevents));
plot(mean_h, c(4), 'kV')
plot(median_h, c(4), 'rV')
xlabel('lags (ms)', 'FontSize', FontSize)
ylabel('Num. of Events', 'FontSize', FontSize)
set(gca, 'FontSize', FontSize)
axis([-100 100 0 c(4)]);
line([0 0], [0 c(4)], 'color', 'k', 'LineWidth', 2, 'LineStyle', '-');
text(.025,.95, strcat('mean= ', num2str(mean_h, '%.1f')), 'color', 'k', 'FontSize', 8, 'Units', 'Normalized')
text(.025,.90, strcat('med.= ', num2str(median_h, '%.1f')), 'color', 'k', 'FontSize', 8, 'Units', 'Normalized')

subplot(3,7,18)
hist(X.max_crosscorr_lag(:,X.Approach_LFPevents), nBins)
c = axis;
hold on
mean_h = mean(X.max_crosscorr_lag(:,X.Approach_LFPevents)); 
median_h = median(X.max_crosscorr_lag(:,X.Approach_LFPevents));
plot(mean_h, c(4), 'kV')
plot(median_h, c(4), 'rV')
xlabel('lags (ms)', 'FontSize', FontSize)
ylabel('Num. of Events', 'FontSize', FontSize)
set(gca, 'FontSize', FontSize)
axis([-100 100 0 c(4)]);
line([0 0], [0 c(4)], 'color', 'k', 'LineWidth', 2, 'LineStyle', '-');
text(.025,.95, strcat('mean= ', num2str(mean_h, '%.1f')), 'color', 'k', 'FontSize', 8, 'Units', 'Normalized')
text(.025,.90, strcat('med.= ', num2str(median_h, '%.1f')), 'color', 'k', 'FontSize', 8, 'Units', 'Normalized')

subplot(3,7,19)
hist(X.max_crosscorr_lag(:,X.Waiting_LFPevents), nBins)
c = axis;
hold on
mean_h = mean(X.max_crosscorr_lag(:,X.Waiting_LFPevents)); 
median_h = median(X.max_crosscorr_lag(:,X.Waiting_LFPevents));
plot(mean_h, c(4), 'kV')
plot(median_h, c(4), 'rV')
xlabel('lags (ms)', 'FontSize', FontSize)
ylabel('Num. of Events', 'FontSize', FontSize)
set(gca, 'FontSize', FontSize)
axis([-100 100 0 c(4)]);
line([0 0], [0 c(4)], 'color', 'k', 'LineWidth', 2, 'LineStyle', '-');
text(.025,.95, strcat('mean= ', num2str(mean_h, '%.1f')), 'color', 'k', 'FontSize', 8, 'Units', 'Normalized')
text(.025,.90, strcat('med.= ', num2str(median_h, '%.1f')), 'color', 'k', 'FontSize', 8, 'Units', 'Normalized')

subplot(3,7,20)
hist(X.max_crosscorr_lag(:,X.Reward_LFPevents), nBins)
c = axis;
hold on
mean_h = mean(X.max_crosscorr_lag(:,X.Reward_LFPevents)); 
median_h = median(X.max_crosscorr_lag(:,X.Reward_LFPevents));
plot(mean_h, c(4), 'kV')
plot(median_h, c(4), 'rV')
xlabel('lags (ms)', 'FontSize', FontSize)
ylabel('Num. of Events', 'FontSize', FontSize)
set(gca, 'FontSize', FontSize)
axis([-100 100 0 c(4)]);
line([0 0], [0 c(4)], 'color', 'k', 'LineWidth', 2, 'LineStyle', '-');
text(.025,.95, strcat('mean= ', num2str(mean_h, '%.1f')), 'color', 'k', 'FontSize', 8, 'Units', 'Normalized')
text(.025,.90, strcat('med.= ', num2str(median_h, '%.1f')), 'color', 'k', 'FontSize', 8, 'Units', 'Normalized')

subplot(3,7,21)
hist(X.max_crosscorr_lag(:,X.PostRun_LFPevents), nBins)
c = axis;
hold on
mean_h = mean(X.max_crosscorr_lag(:,X.PostRun_LFPevents)); 
median_h = median(X.max_crosscorr_lag(:,X.PostRun_LFPevents));
plot(mean_h, c(4), 'kV')
plot(median_h, c(4), 'rV')
xlabel('lags (ms)', 'FontSize', FontSize)
ylabel('Num. of Events', 'FontSize', FontSize)
set(gca, 'FontSize', FontSize)
axis([-100 100 0 c(4)]);
line([0 0], [0 c(4)], 'color', 'k', 'LineWidth', 2, 'LineStyle', '-');
text(.025,.95, strcat('mean= ', num2str(mean_h, '%.1f')), 'color', 'k', 'FontSize', 8, 'Units', 'Normalized')
text(.025,.90, strcat('med.= ', num2str(median_h, '%.1f')), 'color', 'k', 'FontSize', 8, 'Units', 'Normalized')




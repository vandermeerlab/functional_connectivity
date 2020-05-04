function [keep] = copy_of_AUX_shuffle_phases_JS_stored_gamma_times(low_freq, high_freq, varargin)
% shuffles the phases on gamma events that are provided in 3-D matrix X, which is calculated from the function 'find_gamma_bouts.m'
% visualizes the gamma event if doPlot == 1

doPlot = 1;
process_varargin(varargin);

SSN = GetSSN('SingleSession');
if low_freq == 45 && high_freq == 65;
    load(strcat(SSN, '-gamma_bouts-low.mat'));
    
    
    
elseif low_freq == 70 && high_freq == 90;
    load(strcat(SSN, '-gamma_bouts-low.mat'));
else
    error('not sure which gamma events to load')
end


lastEvent = length(gammaIndicesStart);
keep = nan(1,lastEvent);
for iEvent = 1:lastEvent;
    disp(iEvent)
    sig_in = X(1,:,iEvent);
    eeg2 = X(2,:,iEvent); % original vstr time series for gamma event number 1
    md = abs(fft(sig_in));
    phi = angle(fft(sig_in));
    len = floor(length(phi)/2);
    phi_half = phi(1:len);
    %% generate shuffled versions of OFC gamma event and run crosscorrs against VSTR
    numShuffles = 1000;
    max_crosscorr_value_shuff = nan(1,1000);
    for iShuff = 1:numShuffles; % disp(iShuff)
        phi_shuf_half = phi_half(randperm(len)); % permute phases
        if mod(length(phi),2) % odd length
            phi_shuf = cat(2,phi_shuf_half,0,-fliplr(phi_shuf_half));
        else
            phi_shuf = cat(2,phi_shuf_half,-fliplr(phi_shuf_half));
        end
        sig_out = ifft(md.*exp(1i*phi_shuf),'symmetric'); % reconstruct signal with shuffled phases
        eeg1 = sig_out;  % shuffled ofc time series for gamma event number 1
        [~, ~, ~, max_crosscorr_value_shuff(iShuff)]=amp_crosscorr(eeg1, eeg2, 45, 65, samp_freq);
        %     [lags_shuff(iShuff,:), crosscorr_shuff(:,iShuff), max_crosscorr_lag_shuff(iShuff), max_crosscorr_value_shuff(iShuff)]=amp_crosscorr(eeg1, eeg2, 45, 65, samp_freq);
    end
    %% Calculate the REAL crosscorr (nonshuffled) and compare to the shuffled distribution
    eeg1 = sig_in;
    [lags, crosscorr, max_crosscorr_lag, max_crosscorr_value_real, g]=amp_crosscorr(eeg1, eeg2, 45, 65, samp_freq);
    % [lags_real, crosscorr_real, max_crosscorr_lag_real, max_crosscorr_value_real]=amp_crosscorr(eeg1, eeg2, 45, 65, samp_freq);
    percentile = .95;
    sorted_crosscorr_shuff = sort(max_crosscorr_value_shuff);
    percentile_index = percentile*length(sorted_crosscorr_shuff);
    crosscorr_cutoff = sorted_crosscorr_shuff(round(percentile_index));
    if max_crosscorr_value_real >= crosscorr_cutoff;
        keep(iEvent) = 1;
    else
        keep(iEvent) = 0;
    end
    if doPlot == 1;
        clf
        % crosscorr plot fof this gamma event
        subplot(2,2,1)
        plot(lags, crosscorr,'color',[0 0 1],'linewidth',2),hold on %plots crosscorrelations
        plot(lags(g),crosscorr(g),'rp','markerfacecolor',[1 0 0],'markersize',10)%plots marker at the peak of the cross correlation
        plot([0 0],[1.05*max(crosscorr) 0.95*min(crosscorr)],'color',[0 0 0],'linestyle',':', 'linewidth',2) %plots dashed line at zero lag
        set(gca,'xtick',[-100 -50 0 50 100])
        axis tight, box off, xlim([-101 100])
        xlabel('Lag (ms)','fontsize',14)
        ylabel('Crosscorrelation','fontsize',14)
        title(strcat('max crosscorr = ', num2str(max_crosscorr_value_real)))
        % crosscorr values for vstr time series vs. 1,000 phase permuted ofc time series
        subplot(2,2,2)
        hist(max_crosscorr_value_shuff, 50);
        c = axis;
        line([max_crosscorr_value_real max_crosscorr_value_real], [c(3) c(4)], 'Color', 'r', 'LineWidth', 1);
        title(strcat('keep = ', num2str(keep(iEvent))))
        xlabel('max crosscorr value', 'fontsize', 14)
        ylabel('count', 'fontsize', 14)
        % raw LFP traces for OFC and VSTR
        subplot(2,2,3)
        eeg1 = X(1,:,iEvent);
        eeg2 = X(2,:,iEvent);
        hold on
        plot(eeg1, 'r')  % ofc
        plot(eeg2, 'b')  % vstr
        legend('ofc', 'vstr')
        % g50 amplitude envelopes for OFC and VSTR
        subplot(2,2,4)
        Gofc = tsd(lags(1:end-1)', eeg1');
        Gvstr = tsd(lags(1:end-1)', eeg2');
        [~, ~, ~, ~, IEofc]=InstSig(Gofc,45,65,64);
        [~, ~, ~, ~, IEvstr]=InstSig(Gvstr,45,65,256);
        plot(Gofc.range, IEofc.data, 'r')  % ofc
        plot(Gvstr.range, IEvstr.data, 'b')  % vstr
        
        pause
    end
end

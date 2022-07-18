% TR's version of the inverse filtering function
% it can be used for:
% - normalization of HRIR measurements
% - diffuse field equalisation
% - headphone equalisation

function invh_minph = createInverseFilter(h,Fs,Noct,lfefreq)
    Nfft = length(h);
    H = abs(fft(h,Nfft));
    freq = ((0:Nfft-1)*Fs/Nfft)';

    %% spectrum shaping regularization
    if eq(lfefreq,[0 0])
        lfext_amp = 1;
        lfext_smp = 0;
    else
        % establish amplitude level at lf
        idmin = find(freq >= lfefreq(1), 1 );           % f min
        idmax = find(freq <= lfefreq(2), 1, 'last');    % f max
        lfext_amp = 10^(mean(20*log10(H(idmin:idmax)))/20);
        lfext_smp = idmax;
    end

    % establish amplitude level at hf
    hfampmidmin = find(freq >= 16000, 1 );           % f min
    hfampmidmax = find(freq <= 20000, 1, 'last');    % f max
    hfext_amp = 10^(mean(20*log10(H(hfampmidmin:hfampmidmax)))/20);

    % define "crossover" region
    idmin = find(freq >= 18000, 1 );           % f min
    idmax = find(freq <= 20000, 1, 'last');    % f max

    % flatten above fmax
    crossfade_win = hann(2*(idmax-idmin+1));
    win1 = crossfade_win(length(crossfade_win)/2+1:end);
    win2 = crossfade_win(1:length(crossfade_win)/2);
    H_reg = H;
    if lfext_smp ~= 0
        H_reg(1:lfext_smp) = lfext_amp;
    end
    H_reg(idmin:idmax) = H(idmin:idmax) .* win1 + hfext_amp .* win2;
    H_reg(idmax+1:end) = hfext_amp;

    % mirror the left part of H
    H_reg(end/2+1:end) = fliplr(H_reg(2:end/2+1)')';
    
    % freq-domain smoothing using gaussian kernels
    H_sm = smoothSpectrum(20*log10(H_reg(1:end/2)),freq(1:end/2),Noct);
    H_sm = 10.^(H_sm/20);
    H_sm = [H_sm; flipud(H_sm)];
    
    %% calculate inverse filter
    iH=1./H_sm;
    invh=circshift(ifft(iH,'symmetric'),Nfft/2);

    %% kirkeby style regularization and inversion
    freg = [0 18000 20000 freq(end)];
    Lin = -35; % what is below this threshold will be regularised
    Lout = 0; % how much attenuate in dB
    Lreg = [Lin Lin Lout Lout];
    B = interp1(freg,Lreg,freq,'linear');
    B = 10.^(B./20);
    iHkirk=conj(H)./((conj(H).*H)+(conj(B).*B)); % calculating regulated spectral inverse
    
    %% get minimum phase inv h
    invh_minph = minph(invh);
    
    %% PLOTTING
    figure('Name','Inverse Filter','NumberTitle','off','WindowStyle','docked');
    tiledlayout(2,2)
    
%     % impulse responses
%     nexttile
%     hold on
%     box on
%     plot(h)
%     plot(invh-10)
%     plot(invh_minph-20)
%     legend('h','invh','invh_{minph}')

%     % impulse responses
%     nexttile
%     hold on
%     box on
%     sh = length(H)/2;
%     plot(circshift(ifft(H,'symmetric'),sh))
%     plot(circshift(ifft(iH,'symmetric'),sh))
%     plot(circshift(ifft(iHkirk,'symmetric'),sh))
%     yline(0,'--k')
%     xlim([sh sh*1.2]);
%     legend('h','invh_{flattening}','invh_{kirkeby}')

    % magnitudes
    nexttile
    hold on
    box on
    plot(freq,20*log10(H))
%     plot(freq,20*log10(H_reg))
%     plot(freq,20*log10(H_sm))
    plot(freq,20*log10(iH))
    plot(freq,20*log10(iHkirk))
    plot(freq,20*log10(B),'--k')
%     legend('H','H_{reg}','H_{sm}','iH','iH_{kirk}','B')
    legend('Measured magnitude','Inversed using spectrum shaping','Inversed using Kirkeby regularization','Kirkeby regularization parameter')
    set(gca,'xscale','log')
    xlim([20 Fs/2]);
    ylim([-40 40]);
    xlabel('Frequency (Hz)')
    ylabel('Magnitude (dB)')

    % inversed magnitudes (HF)
    nexttile
    hold on
    box on
    plot(freq,20*log10(H))
    plot(freq,20*log10(iH))
    plot(freq,20*log10(iHkirk))
    legend('Measured magnitude','Inversed using spectrum shaping','Inversed using Kirkeby regularization')
    set(gca,'xscale','log')
    xlim([12000 Fs/2]);
    ylim([-40 40]);
    xlabel('Frequency (Hz)')
    ylabel('Magnitude (dB)')

    % magnitudes
    nexttile
    hold on
    box on
    plot(freq,20*log10(H))
    plot([freq(hfampmidmin) freq(hfampmidmax)],[20*log10(hfext_amp) 20*log10(hfext_amp)],'--','LineWidth',2)
    plot(freq(idmin:idmax),20*log10(H(idmin:idmax) .* win1))
    plot(freq(idmin:idmax),20*log10(hfext_amp .* win2))
    plot(freq,20*log10(H_reg),'LineWidth',2)
    xline(freq(idmin),'--k')
    xline(freq(idmax),'--k')
    legend('Measured magnitude','Mean level at HF','Truncated magnitude','HF extension','Regularized magnitude','Location','northwest')
%     text(18060,10,'Cross-fade region')
    set(gca,'xscale','log')
    xlim([12000 Fs/2]);
    ylim([-40 40]);
    xlabel('Frequency (Hz)')
    ylabel('Magnitude (dB)')

    % save figure
    figlen = 12;
    width = 3*figlen;
    height = 2*figlen;
    set(gcf,'Units','centimeters','PaperPosition',[0 0 width height],'PaperSize',[width height]);
    saveas(gcf,'figures/normalization_flattening.pdf')
%     close


end
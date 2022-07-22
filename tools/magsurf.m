function magsurf(Fs,ang,hrir,ytit,tit,zlimit)
    Nfft = 1024;
    f = ((0:Nfft-1)*Fs/Nfft)';
    magnitude = abs(fft(hrir,Nfft));
    magnitude = 20*log10(magnitude);
    s = pcolor(f,ang,magnitude');
    s.EdgeColor = 'none';
    hold on
    yline(-90,'--')
    yline(0,'--')
    yline(90,'--')
    set(gca,'xscale','log')
    xlim([500 24000])
    xticks([1000 2000 4000 8000 16000])
    ylim([-180 180])
    zlim(zlimit)
    caxis(zlimit)
    xlabel('Frequency (Hz)')
    ylabel(ytit)
    title(tit)
    colorbar
end

% function magsurf(Fs,ang,hrir,ytit,tit,zlimit)
%     Nfft = 1024;
%     f = ((0:Nfft-1)*Fs/Nfft)';
%     magnitude = abs(fft(hrir,Nfft));
%     magnitude = 20*log10(magnitude);
%     s = pcolor(f,ang,magnitude');
%     s.EdgeColor = 'none';
%     hold on
%     yline(-90,'--')
%     yline(0,'--')
%     yline(90,'--')
%     xlim([1000 12000])
%     xticks([1000 2000 4000 8000 16000])
%     ylim([-180 180])
%     zlim(zlimit)
%     caxis(zlimit)
%     xlabel('Frequency (Hz)')
%     ylabel(ytit)
%     title(tit)
%     colorbar
% end
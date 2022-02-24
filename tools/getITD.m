% get ITD, direct sound sample indices, and direct sound delay

function [ITD, maxL, maxR, dlyL, dlyR] = getITD(irLeft, irRight, Fs)

    fx = 1500/(Fs/2);
    f = [0 fx fx 1];
    mhi = [1 1 0 0];
    B = fir2(8,f,mhi);
    A = 1;
    irLeft = filter(B, A, irLeft);
    irRight = filter(B, A, irRight);
    gd = grpdelay(B,A);
    filter_dly = median(gd(1:end/2));
%     filter_dly = 0;

    % find peak value and use windowing to cut off strong room
    % reflections in the contralateral ear
    [~,miridx_left] = max(abs(irLeft));
    [~,miridx_right] = max(abs(irRight));
    
    peakidx = min(miridx_left,miridx_right);
    winstart = peakidx-(Fs*0.0005);
    winend = peakidx+(Fs*0.002);
    
%     if true
%         figure('Name','ITD debug','NumberTitle','off','WindowStyle','docked')
%         subplot(2,1,1)
%         hold on
%         plot(irLeft)
%         xline(peakidx,'--r')
%         xline(winstart,'--k')
%         xline(winend,'--k')
%         xlim([peakidx - 300 peakidx + 300])
% 
%         subplot(2,1,2)
%         hold on
%         plot(irRight)
%         xline(peakidx,'--r')
%         xline(winstart,'--k')
%         xline(winend,'--k')
%         xlim([peakidx - 300 peakidx + 300])
%     end

    % upsampling
    r = 8;
    irLeftInterp = interp(irLeft(winstart:winend),r);
    irRightInterp = interp(irRight(winstart:winend),r);        

    % left ear delay
    [acor, lag] = xcorr(irLeftInterp,minph(irLeftInterp));
    [~,indexL] = max(acor);
    dly_smp_left = lag(indexL);

    % right ear delay
    [acor, lag] = xcorr(irRightInterp,minph(irRightInterp));
    [~,indexR] = max(acor);
    dly_smp_right = lag(indexR);

    ITD = (dly_smp_right - dly_smp_left) * 10^6 / (Fs*r);
    maxL = winstart-1 + round(dly_smp_left / r) - filter_dly;
    maxR = winstart-1 + round(dly_smp_right / r) - filter_dly;
    dlyL = (dly_smp_left / r + (winstart-1 - filter_dly)) * 10^6 / Fs;
    dlyR = (dly_smp_right / r + (winstart-1 - filter_dly)) * 10^6 / Fs;

    if ITD > 1000
        figure('Name','ITD debug','NumberTitle','off','WindowStyle','docked')
        subplot(2,1,1)
        hold on
        plot(irLeft)
        plot(minph(irLeft))
        xline(peakidx,'--r')
        xline(winstart,'--k')
        xline(winend,'--k')

        subplot(2,1,2)
        hold on
        plot(irRight)
        plot(minph(irRight))
        xline(peakidx,'--r')
        xline(winstart,'--k')
        xline(winend,'--k')
    end
 
%     function [h_min] = minph(h)
%         % based on the invFIR function - look there for explanations
%         % https://uk.mathworks.com/matlabcentral/fileexchange/19294-inverse-fir-filter
%         n = length(h);
%         h_cep = real(ifft(log(abs(fft(h(:,1))))));
%         odd = fix(rem(n,2));
%         wn = [1; 2*ones((n+odd)/2-1,1) ; ones(1-rem(n,2),1); zeros((n+odd)/2-1,1)];
%         h_min = zeros(size(h(:,1)));
%         h_min(:) = real(ifft(exp(fft(wn.*h_cep(:)))));
%     end
end
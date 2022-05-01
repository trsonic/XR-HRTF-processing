% get ITD, direct sound sample indices, and direct sound delay

function [ITD, maxL, maxR, dlyL, dlyR] = getITD(irLeft, irRight, Fs)
    fx = 1500/(Fs/2);
    f = [0 fx fx 1];
    mhi = [1 1 0 0];
    B = fir2(8,f,mhi);
    A = 1;
    irLeftFiltered = filter(B, A, irLeft);
    irRightFiltered = filter(B, A, irRight);

    % find peak value and use windowing to cut off strong room
    % reflections in the contralateral ear
    [~,miridx_left] = max(abs(irLeftFiltered));
    [~,miridx_right] = max(abs(irRightFiltered));
    
    peakidx = min(miridx_left,miridx_right);
    winstart = peakidx-fix(Fs*0.001);
    winend = peakidx+fix(Fs*0.003);

    % upsampling
    r = 8;
    irLeftInterp = interp(irLeftFiltered(winstart:winend),r);
    irRightInterp = interp(irRightFiltered(winstart:winend),r);        

    % left ear delay
    [acor, lag] = xcorr(irLeftInterp,minph(irLeftInterp));
    [~,indexL] = max(acor);
    dly_smp_left = lag(indexL);

    % right ear delay
    [acor, lag] = xcorr(irRightInterp,minph(irRightInterp));
    [~,indexR] = max(acor);
    dly_smp_right = lag(indexR);

    ITD = (dly_smp_right - dly_smp_left) * 10^6 / (Fs*r);
    maxL = winstart-1 + round(dly_smp_left / r);
    maxR = winstart-1 + round(dly_smp_right / r);
    dlyL = (dly_smp_left / r + (winstart-1)) * 10^6 / Fs;
    dlyR = (dly_smp_right / r + (winstart-1)) * 10^6 / Fs;

    if ITD > 1000
        figure('Name','ITD debug','NumberTitle','off','WindowStyle','docked')
        subplot(2,1,1)
        hold on
        plot(irLeft)
        plot(irLeftFiltered)
        xline(peakidx,'--r')
        xline(winstart,'--k')
        xline(winend,'--k')
        xline(maxL,'-g')
        xlim([peakidx-200 peakidx+600])

        subplot(2,1,2)
        hold on
        plot(irRight)
        plot(irRightFiltered)
        xline(peakidx,'--r')
        xline(winstart,'--k')
        xline(winend,'--k')
        xline(maxR,'-g')
        xlim([peakidx-200 peakidx+600])
    end
end
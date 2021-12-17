function [ITD, maxL, maxR] = getITD(irLeft, irRight, Fs)

    % filter
    xfreq = 1500;
    filter_order = 2;    
    [B,  A ] = butter( filter_order, xfreq/Fs*2, 'low'  );
    irLeft = filter(B, A, irLeft);
    irRight = filter(B, A, irRight);
    irLeft = filter(B, A, irLeft);
    irRight = filter(B, A, irRight);

    % find peak value and use windowing to cut off strong room
    % reflections in the contralateral ear
    [~,miridx_left] = max(abs(irLeft));
    [~,miridx_right] = max(abs(irRight));

    peakidx = min(miridx_left,miridx_right);
    winstart = peakidx-(Fs*0.0005);
    winend = peakidx+(Fs*0.003);

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
    maxL = winstart-1 + round(dly_smp_left / r);
    maxR = winstart-1 + round(dly_smp_right / r);

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
 
    function [h_min] = minph(h)
        % based on the invFIR function - look there for explanations
        % https://uk.mathworks.com/matlabcentral/fileexchange/19294-inverse-fir-filter
        n = length(h);
        h_cep = real(ifft(log(abs(fft(h(:,1))))));
        odd = fix(rem(n,2));
        wn = [1; 2*ones((n+odd)/2-1,1) ; ones(1-rem(n,2),1); zeros((n+odd)/2-1,1)];
        h_min = zeros(size(h(:,1)));
        h_min(:) = real(ifft(exp(fft(wn.*h_cep(:)))));
    end
end
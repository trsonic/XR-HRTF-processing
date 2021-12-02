function [ITD, maxL, maxR] = getITD(irLeft, irRight, Fs)
        % find peak value and use windowing to cut off strong room
        % reflections in the contralateral ear
        [~,miridx_left] = max(abs(irLeft));
        [~,miridx_right] = max(abs(irRight));

        peakidx = min(miridx_left,miridx_right);
        winstart = peakidx-(Fs*0.0005);
        winend = peakidx+(Fs*0.001);
       
        % upsampling
        r = 8;
        irLeftInterp = interp(irLeft(winstart:winend),r);
        irRightInterp = interp(irRight(winstart:winend),r);

        % left ear delay
        [acor, lag] = xcorr(irLeftInterp,minph(irLeftInterp));
%         [~,index] = max(abs(acor));
        [~,indexL] = max(acor);
        dly_smp_left = lag(indexL);

        % right ear delay
        [acor, lag] = xcorr(irRightInterp,minph(irRightInterp));
%         [~,index] = max(abs(acor));
        [~,indexR] = max(acor);
        dly_smp_right = lag(indexR);

        ITD = (dly_smp_right - dly_smp_left) * 10^6 / (Fs*r);
        maxL = winstart-1 + round(dly_smp_left / r);
        maxR = winstart-1 + round(dly_smp_right / r);

        if ITD > 1000
            figure('Name','ITD debug','NumberTitle','off','WindowStyle','docked')
            hold on
            plot(irLeft)
            plot(minph(irLeft))
            plot(irRight)
            plot(minph(irRight))
            xline(peakidx,'--k')
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
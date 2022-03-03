% TR's version of the inverse filtering function
function invh_minph = createInverseFilter(h,Fs,Noct,lfefreq)

%     if exist('input_type', 'var')
%         if strcmp(input_type,'h')
%             Nfft = length(h);
%             H=abs(fft(h,Nfft));
%         elseif strcmp(input_type,'H')
%             Nfft = length(h);
%             H=h;
%         end
%     else
%         Nfft = length(h);
%         H=abs(fft(h,Nfft));
%     end
    
    Nfft = length(h);
    H=abs(fft(h,Nfft));

    disp(['length: ' num2str(length(h)) ', Fs: ' num2str(Fs)])
    
    %% regularization
    H_reg = regularizeH(H,Fs, lfefreq);
    
    %% freq-domain smoothing if necessary...
    freq = ((0:Nfft-1)*Fs/Nfft)';
    H_sm = smoothSpectrum(H_reg(1:end/2),freq(1:end/2),Noct);
    H_sm = [H_sm; fliplr(H_sm')'];
    
    %% calculate inverse filter
    iH=conj(H_sm)./(conj(H_sm).*H_sm);
%     iH=conj(H_reg)./(conj(H_reg).*H_reg);
%     iH=conj(H)./((conj(H).*H)+(conj(B).*B)); % calculating regulated spectral inverse
    invh=circshift(ifft(iH,'symmetric'),Nfft/2);
    
%     %% window the inversed h
%     L = Nfft;
%     win = hann(L); %0.5*(1-cos(2*pi*(1:L)'/(L+1)));
%     invh_win=win.*invh(end/2-L/2+1:end/2+L/2); % truncation to length L
    
    %% get minimum phase inv h
    invh_minph = minph(invh);
    
%     %% PLOTTING
%     figure('Name','Inverse Filter','NumberTitle','off','WindowStyle','docked');
%     tiledlayout(1,2)
%     nexttile
%     hold on
%     box on
%     plot(h)
%     plot(invh-10)
%     plot(invh_minph-20)
% 
%     nexttile
%     hold on
%     box on
%     plot(20*log10(H))
%     plot(20*log10(H_reg))
%     plot(20*log10(iH))
    


    function H_fixed = regularizeH(H,Fs, lfefreq)
        f = (0:length(H)-1)*Fs/length(H);
        
        if ~eq(lfefreq,[0 0])
            % establish amplitude level at lf
            idmin = find(f >= lfefreq(1), 1 );           % f min
            idmax = find(f <= lfefreq(2), 1, 'last');    % f max
    %         idmin = find(f >= 60, 1 );           % f min
    %         idmax = find(f <= 300, 1, 'last');    % f max
            lfext_amp = 10^(mean(20*log10(H(idmin:idmax)))/20);
            lfext_smp = idmax;
        else
            lfext_amp = 1;
            lfext_smp = 0;
        end

        % establish amplitude level at hf
        idmin = find(f >= 16000, 1 );           % f min
        idmax = find(f <= 20000, 1, 'last');    % f max
        hfext_amp = 10^(mean(20*log10(H(idmin:idmax)))/20);


        % define "crossover" region
        idmin = find(f >= 18000, 1 );           % f min
        idmax = find(f <= 20000, 1, 'last');    % f max

        % flatten above fmax
        crossfade_win = hann(2*(idmax-idmin+1));
        win1 = crossfade_win(length(crossfade_win)/2+1:end);
        win2 = crossfade_win(1:length(crossfade_win)/2);
        H_fixed = H;
        if lfext_smp ~= 0
            H_fixed(1:lfext_smp) = lfext_amp;
        end
        H_fixed(idmin:idmax) = H(idmin:idmax) .* win1 + hfext_amp .* win2;
        H_fixed(idmax+1:end) = hfext_amp;

        % mirror the left part of H
        H_fixed(end/2+1:end) = fliplr(H_fixed(2:end/2+1)')';
    end
end



% function ir_fixed = flattenMagHF(ir, Fs)
%     Nfft = length(ir);
%     f = (Fs/Nfft:Fs/Nfft:Fs);
%     mag_target = abs(fft(ir, Nfft));
%     
%     % establish amplitude level
%     idmin = find(f >= 16000, 1 );           % f min
%     idmax = find(f <= 20000, 1, 'last');    % f max
%     hfext_amp = mean(mag_target(idmin:idmax));
%     
%     % define "crossover" region
%     idmin = find(f >= 18000, 1 );           % f min
%     idmax = find(f <= 20000, 1, 'last');    % f max
%     
%     % flatten above fmax
%     win = hann(2*(idmax-idmin+1));
%     win1 = win(length(win)/2+1:end);
%     win2 = win(1:length(win)/2);
%     mag_fixed = mag_target;
%     mag_fixed(idmin:idmax) = mag_target(idmin:idmax) .* win1 + hfext_amp .* win2;
%     mag_fixed(idmax+1:end) = hfext_amp;
%     
%     % back to time domain
%     ir_fixed = ifft(mag_fixed,'symmetric');
%     ir_fixed=circshift(ir_fixed,Nfft/2,1);
% %     ir_fixed=0.5*(1-cos(2*pi*(1:Nfft)'/(Nfft+1))).*ir_fixed;
% end
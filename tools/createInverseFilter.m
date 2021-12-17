% TR's version of the inverse filtering function
function invh_minph = createInverseFilter(h,Fs)
    disp(['length: ' num2str(length(h)) ', Fs: ' num2str(Fs)])
    
    Nfft = length(h);
    H=abs(fft(h,Nfft));
    
    %% regularization
    H_reg = flattenMagHF(H,Fs);
    
    %% freq-domain smoothing if necessary...
    
    %% calculate inverse filter
    iH=conj(H_reg)./(conj(H_reg).*H_reg);
%     iH=conj(H)./((conj(H).*H)+(conj(B).*B)); % calculating regulated spectral inverse
    invh=circshift(ifft(iH,'symmetric'),Nfft/2);
    
%     %% window the inversed h
%     L = Nfft;
%     win = hann(L); %0.5*(1-cos(2*pi*(1:L)'/(L+1)));
%     invh_win=win.*invh(end/2-L/2+1:end/2+L/2); % truncation to length L
    
    %% get minimum phase inv h
    invh_minph = minph(invh);
    
    %% PLOTTING
    figure('Name','Inverse Filter','NumberTitle','off','WindowStyle','docked');
    tiledlayout(1,2)
    nexttile
    hold on
    box on
    plot(h)
    plot(invh-10)
    plot(invh_minph-20)

    nexttile
    hold on
    box on
    plot(20*log10(H))
    plot(20*log10(H_reg))
    plot(20*log10(iH))
    
    % calculate minimum phase component of impulse response
    function [h_min] = minph(h)
        n = length(h);
        h_cep = real(ifft(log(abs(fft(h(:,1))))));
        odd = fix(rem(n,2));
        wn = [1; 2*ones((n+odd)/2-1,1) ; ones(1-rem(n,2),1); zeros((n+odd)/2-1,1)];
        h_min = zeros(size(h(:,1)));
        h_min(:) = real(ifft(exp(fft(wn.*h_cep(:)))));
    end

    function H_fixed = flattenMagHF(H,Fs)
        f = (0:length(H)-1)*Fs/length(H);
        
        % establish amplitude level at hf
        idmin = find(f >= 60, 1 );           % f min
        idmax = find(f <= 300, 1, 'last');    % f max
        lfext_amp = 10^(mean(20*log10(H(idmin:idmax)))/20);
        lfext_smp = idmax;

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
        H_fixed(1:lfext_smp) = lfext_amp;
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
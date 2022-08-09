function hpEQ(hpirBank, subjectdir)
    % calculate average magnitude for left & right   
    for i = 1:length(hpirBank)
        if i == 1
            Nfft = size(hpirBank(i).fullIR,1);
            mag_ir_avgL = zeros(Nfft,1);
            mag_ir_avgR = zeros(Nfft,1);
            sn = 1 / length(hpirBank);
        end
        mag_ir_avgL = mag_ir_avgL + (abs(fft(hpirBank(i).fullIR(:,1), Nfft)).^2) * sn;
        mag_ir_avgR = mag_ir_avgR + (abs(fft(hpirBank(i).fullIR(:,2), Nfft)).^2) * sn;
    end
    
    mag_ir_avgL = sqrt(mag_ir_avgL);
    mag_ir_avgR = sqrt(mag_ir_avgR);

    % back to time domain
    ir_avgL = ifft(mag_ir_avgL,'symmetric');
    ir_avgR = ifft(mag_ir_avgR,'symmetric');
    ir_avgL = circshift(ir_avgL,Nfft/2,1);
    ir_avgR = circshift(ir_avgR,Nfft/2,1);
    hpir = [ir_avgL ir_avgR];
    Fs = hpirBank(1).Fs;

    % adjust levels
    [f,mag] = getMagnitude(hpir,Fs,'log');
    idmin = find(f >= 200, 1 );           % f min
    idmax = find(f <= 800, 1, 'last');    % f max
    gain1 = mean(mag(idmin:idmax,1));
    gain2 = mean(mag(idmin:idmax,2));
    hpir(:,1) = hpir(:,1) * 10^(-gain1/20);
    hpir(:,2) = hpir(:,2) * 10^(-gain2/20);
    
    % create EQ filters
    invh(:,1) = createInverseFilter(hpir(:,1), Fs, 12, [400 600]);
    invh(:,2) = createInverseFilter(hpir(:,2), Fs, 12, [400 600]);

    %% plotting
    figure('Name','HPTF + inv resp','NumberTitle','off','WindowStyle','docked');
    hold on
    box on

    [f,mag] = getMagnitude(hpir(:,1),Fs,'log');
    plot(f,mag,'-g','LineWidth',1);

    [f,mag] = getMagnitude(hpir(:,2),Fs,'log');
    plot(f,mag,'-r','LineWidth',1);

    [f,mag] = getMagnitude(invh(:,1),Fs,'log');
    plot(f,mag,'--g','LineWidth',2);

    [f,mag] = getMagnitude(invh(:,2),Fs,'log');
    plot(f,mag,'--r','LineWidth',2);

    set(gca,'xscale','log')
    grid on
    xlim([20 Fs/2]);
    ylim([-20 20]);
    %             legend('Left channel', 'Right channel','','', 'Left inverse filter', 'Right inverse filter','location','northwest')
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    
    figure('Name','HPIR HPEQIR','NumberTitle','off','WindowStyle','docked');
    hold on
    plot(hpir(:,1))
    plot(hpir(:,2))
    plot(invh(:,1))
    plot(invh(:,2))
    
    mkdir([subjectdir '/hpeq/'])
    audiowrite([subjectdir '/hpeq/' 'hpeq.wav'], invh * 0.25, Fs)
end
close all
clear

addpath("tools/")

load("data\20220807-XR-SUBJ001B\irBankProcessed.mat")


figure('Name','HRTF magnitude','NumberTitle','off','WindowStyle','docked');

for i = 3 %1:length(irBank)
    Fs = irBank(i).Fs;
    hold off
    [f,mag] = getMagnitude(irBank(i).rawHRIR(:,1),Fs,'log');
    plot(f,mag,'-b','LineWidth',1);
    hold on
    [f,mag] = getMagnitude(irBank(i).rawHRIR(:,2),Fs,'log');
    plot(f,mag,'-r','LineWidth',1);
    

    legend('Left Ear','Right Ear','location','southwest')
    xticks([125 250 500 1000 2000 4000 8000 16000])
    xlim([200 18000]);
    ylim([-30 30]);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    title(['Azimuth: ' num2str(irBank(i).azimuth) '^{\circ}, elevation: ' num2str(irBank(i).elevation) '^{\circ}, distance: ' num2str(irBank(i).distance) ' m']);
    set(gca,'xscale','log')

%     pause

end
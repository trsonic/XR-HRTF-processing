close all
clear


%% kirkeby B parameter
figure
Nfft = 512;
Fs = 48000;
freq = ((0:Nfft-1)*Fs/Nfft)';

freg = [0 160 250 16000 18000 Fs];
Lreg = [120 120 -160 -160 120 120];
B = interp1(freg,Lreg,freq,"linear");
B = 10.^(-B./20);

plot(freq,B)




%% hanning window crossfade
figure
hold on

crossfade_win = hann(512);
win1 = crossfade_win(length(crossfade_win)/2+1:end);
win2 = crossfade_win(1:length(crossfade_win)/2);

win2 = win2*0.5;

% plot(crossfade_win)
plot(win1)
plot(win2)

plot(win1+win2)

ylim([-0.1 1.1])




figure

subplot(2,1,1)
hold on
ylim([0 1])
ylabel('Gain (linear)')

for p = 1:6
    win = getWindow(p);
    plot(win)
end

subplot(2,1,2)
hold on
ylim([-60 0])
ylabel('Gain (dB)')
for p = 1:6
    win = getWindow(p);
    plot(20*log10(win))
end



function win = getWindow(p)
    % define window    
    win1 = hann(80).^p;
    win1 = win1(1:end/2);
    win2 = hann(200).^p;
    win2 = win2(end/2+1:end);
    win = [win1; ones(40,1); win2;];
    winshift = length(win1); % how many samples the window should be shifted forward from the peak

end
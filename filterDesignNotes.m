close all
clear


    % filter
%     xfreq = 1500;
%     filter_order = 2;    
%     [B,  A ] = butter( filter_order, xfreq/Fs*2, 'low'  );
%     irLeft = filter(B, A, irLeft);
%     irRight = filter(B, A, irRight);


%    lpFilt = designfilt('lowpassfir', 'PassbandFrequency', 1200,...
%            'StopbandFrequency', 1500, 'PassbandRipple', 0.5, ...
%            'StopbandAttenuation', 65, 'DesignMethod', 'kaiserwin','SampleRate',Fs);
                     
    
%     B = fir1(34,0.0625,chebwin(35,30));
Fs = 48000;
fx = 1500/(Fs/2);
f = [0 fx fx 1];
mhi = [1 1 0 0];
B = fir2(8,f,mhi);
A = 1;

freqz(B,A,2048,Fs)

h = firminphase(B);

grpdelay(B,A,2048)
    
%     Fs = 48000;
%     disp(median(grpdelay(B,A)));

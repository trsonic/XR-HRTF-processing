function [f,mag] = getMagnitude(ir,Fs,linlog)
    Nfft = length(ir);
    f = (Fs/Nfft:Fs/Nfft:Fs)';
   
    mag = abs(fft(ir, Nfft));
    
    if strcmp(linlog,'log')
        mag = 20*log10(mag);
    end
end

function [f,mag] = getMagnitude(ir,Fs,linlog)
    Nfft = length(ir);
    f = ((0:Nfft-1)*Fs/Nfft)';
    mag = abs(fft(ir, Nfft));
    
    if strcmp(linlog,'log')
        mag = 20*log10(mag);
    end
end


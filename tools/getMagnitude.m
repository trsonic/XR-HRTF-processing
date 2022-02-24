function [f,mag] = getMagnitude(ir,Fs,linlog)
for i = 1:size(ir,2)
    Nfft = size(ir,1);
    f = ((0:Nfft-1)*Fs/Nfft)';
    mag(:,i) = abs(fft(ir(:,i), Nfft));
    
    if strcmp(linlog,'log')
        mag(:,i) = 20*log10(mag(:,i));
    end
end
end


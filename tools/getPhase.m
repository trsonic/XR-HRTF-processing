function [f,ph] = getPhase(ir,Fs)
    for i = 1:size(ir,2)
        Nfft = size(ir,1);
        f = ((0:Nfft-1)*Fs/Nfft)';
        ph(:,i) = angle(fft(ir(:,i), Nfft));
    end
end
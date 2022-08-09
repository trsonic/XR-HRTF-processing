function ir = deconvolve(sweep,inv_sweep)
    N1 = size(sweep,1) + size(inv_sweep,1);
    ir(:,1) = ifft((fft(sweep(:,1),N1)).*fft(inv_sweep,N1));
    
    if size(sweep,2) == 2
        ir(:,2) = ifft((fft(sweep(:,2),N1)).*fft(inv_sweep,N1));
    end
end
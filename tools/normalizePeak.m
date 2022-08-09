function irBank = normalizePeak(irBank)
    max_val = [];
    
    for i = 1:length(irBank)
        max_val = [max_val; max(max(abs(irBank(i).fullIR)))];
    end
    
    max_val = max(max_val);

    for i = 1:length(irBank)
        irBank(i).fullIR = irBank(i).fullIR ./ max_val;

    end
end
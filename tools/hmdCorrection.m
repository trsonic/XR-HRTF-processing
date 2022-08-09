function irBank = hmdCorrection(irBank)
    load('data/model_interp.mat')
    
    for i = 1:length(irBank)
        if irBank(i).ref == 0
            dist = distance(irBank(i).elevation,irBank(i).azimuth,[model_interp.el],[model_interp.az]);
            [~,idxl] = min(dist);
            dist = distance(irBank(i).elevation,irBank(i).azimuth,[model_interp.el],-[model_interp.az]);
            [~,idxr] = min(dist);
            
            % correct magnitude
            left = irBank(i).fullIR(:,1);
            right = irBank(i).fullIR(:,2);
            left = conv(model_interp(idxl).invh,left);
            right = conv(model_interp(idxr).invh,right); 

            % correct time of arrival
            dly = model_interp(idxl).dtoa_diff;
            shift = -1 * dly * 10^-6 * irBank(i).Fs;
            left = fraccircshift(left,shift);
            
            dly = model_interp(idxr).dtoa_diff;
            shift = -1 * dly * 10^-6 * irBank(i).Fs;
            right = fraccircshift(right,shift);
                  
            irBank(i).fullIR = [left right];    
        end
    end
end
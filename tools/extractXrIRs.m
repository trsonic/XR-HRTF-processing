function irBank = extractXrIRs(subjectdir,inv_pol)
    sweepdir = [subjectdir 'sweeps/'];
    
    %% HRIRS + reference
    [y_inv_sweep, Fs] = audioread([sweepdir 'LS_inv_sweep.wav']);
    filelist = dir([sweepdir '*_measured_*.wav']);
    filelist = [filelist; dir([sweepdir '00_reference_01.wav'])];
    
    irBank = struct;
    for i = 1:length(filelist)
        [irBank(i).y, irBank(i).Fs] = audioread([filelist(i).folder '/' filelist(i).name]);
        irBank(i).name = filelist(i).name;

        if inv_pol
            irBank(i).y = -irBank(i).y;
        end
        name = extractAfter(irBank(i).name,'target');
        irBank(i).azimuth = str2double(extractBefore(extractAfter(name,'_az_'),'_el_'));
        irBank(i).elevation = str2double(extractBefore(extractAfter(name,'_el_'),'_dist_'));
        name = extractAfter(irBank(i).name,'measured');
        irBank(i).distance = str2double(extractBefore(extractAfter(name,'_dist_'),'.wav'));

        if contains(irBank(i).name,'reference','IgnoreCase',true)
            irBank(i).ref = true;
        else
            irBank(i).ref = false;
        end
        irBank(i).lspk = 1;
        
        if irBank(i).Fs ~= Fs
            disp('Fs mismatch!!!')
        end
        irBank(i).fullIR = deconvolve(irBank(i).y,y_inv_sweep);
        
        % skip the first half of the ir
        cut_start = floor(size(irBank(i).fullIR,1)/2);
        cut_end = cut_start + 8191;
        irBank(i).fullIR = irBank(i).fullIR(cut_start:cut_end,:);
    end
    irBank = rmfield(irBank, 'y');
    
%     % HpTF
%     [y_inv_sweep, Fs] = audioread([sweepdir 'HP_inv_sweep.wav']);
%     filelist = dir([sweepdir '00_hptf_*.wav']);
%     
%     hpirBank = struct;
%     for i = 1:length(filelist)
%         [hpirBank(i).y, hpirBank(i).Fs] = audioread([filelist(i).folder '/' filelist(i).name]);
%         hpirBank(i).name = filelist(i).name;
%         if hpirBank(i).Fs ~= Fs
%             disp('Fs mismatch!!!')
%         end
%         hpirBank(i).fullIR = deconvolve(hpirBank(i).y,y_inv_sweep);
%         
%         % skip the first half of the ir
%         cut_start = floor(size(hpirBank(i).fullIR,1)/2);
%         cut_end = cut_start + 8191;
%         hpirBank(i).fullIR = hpirBank(i).fullIR(cut_start:cut_end,:);
%     end
%     hpirBank = rmfield(hpirBank, 'y');
    
end
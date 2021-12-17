close all
clear

% %% RIG
% dirlist = dir('data/*-RIG-*');
% for i = 1:length(dirlist)
%     subjectdir = [dirlist(i).folder '/' dirlist(i).name '/']; inv_pol = false;
%     irBank = extractRigIRs(subjectdir,inv_pol,false);
%     irBank = adjustGains(irBank,-18,-30);
%     irBank = normalizePeak(irBank);
%     plotAndSave(subjectdir,irBank);
% end

%% XR
% subjectdir = 'data/20201217-122pt-2.5m-dayton_vt/'; inv_pol = true;
% subjectdir = 'data/20201217-122pt-2.5m-canford_vt/'; inv_pol = true;
% subjectdir = 'data/20211012-q2_tr/'; inv_pol = true;
% subjectdir = 'data/20211105-A-Jan/'; inv_pol = true;
% subjectdir = 'data/20211126-XR-TR/'; inv_pol = false;
subjectdir = 'data/20211126-XR-Gavin/'; inv_pol = false;
irBank = extractXrIRs(subjectdir,inv_pol);
irBank = normalizePeak(irBank);
plotAndSave(subjectdir,irBank);

function irBank = extractRigIRs(subjectdir,inv_pol,plotting)
    sweepdir = [subjectdir 'sweeps/'];
    [y_inv_sweep, Fs] = audioread([sweepdir 'ZZ_inv_sweep.wav']);

    filelist = dir([sweepdir '*-*.wav']);

    irBank = struct;
    for i = 1:length(filelist)
        disp(i)
        [y, Fs] = audioread([filelist(i).folder '/' filelist(i).name]);
        silence_length = 96000;
        msweep_length = 1999200;
        y=y(silence_length+1:silence_length+msweep_length,:);
        if inv_pol
            y = -y;
        end
        
        irArray = deconvolve(y,y_inv_sweep);
        
        if plotting
            figure('Name','Measured IRs','NumberTitle','off','WindowStyle','docked');
            hold on
            plot(irArray)
        end
        
        msweep_offset = floor(0.6*48000);
        start_offset = floor(length(irArray)/4);
        num_lspk = 50;
        for lspk = 1:num_lspk
            ir_center = start_offset + msweep_offset * (lspk-3);
            cut_start = ir_center - 2048;
            cut_end = cut_start + 8191;
            if plotting
                xline(ir_center, '--')
                xline(cut_start, 'g--')
                xline(cut_end, 'g--')
            end
            
            % write the info
            j = (i-1)*num_lspk + lspk; 
            deg = str2double(extractAfter(extractBefore(filelist(i).name,'deg'),'KEMAR-'));
            [azimuth, elevation] = getLspkAzel(lspk);
            if isnan(deg)
                irBank(j).azimuth = azimuth;
            else
                irBank(j).azimuth = azimuth - deg;
            end
            irBank(j).elevation = elevation;
            irBank(j).distance = 1.5;
            irBank(j).Fs = Fs;
            irBank(j).name = filelist(i).name;
            if contains(filelist(i).name,'reference','IgnoreCase',true)
                irBank(j).ref = true;
            else
                irBank(j).ref = false;
            end
            
            irBank(j).lspk = lspk;
            % write single HRIRs
            irBank(j).fullIR = irArray(cut_start:cut_end,:);
            
            % duplicate left channel to right channel if it's missing (for
            % the reference measurement only)
            if irBank(j).ref && size(irBank(j).fullIR,2) == 1
                irBank(j).fullIR(:,2) = irBank(j).fullIR;
            end
        end
    end
    
    function [az, el] = getLspkAzel(lspk)
        % azel_rig = [0,90;45,65;135,65;-135,65;-45,65;0,45;90,45;180,45;-90,45;45,35;135,35;-135,35;-45,35;18,18;72,18;108,18;162,18;-162,18;-108,18;-72,18;-18,18;0,0;45,0;90,0;135,0;180,0;-135,0;-90,0;-45,0;18,-18;72,-18;108,-18;162,-18;-162,-18;-108,-18;-72,-18;-18,-18;45,-35;135,-35;-135,-35;-45,-35;0,-45;90,-45;180,-45;-90,-45;45,-65;135,-65;-135,-65;-45,-65;0,-90];
        azel_rig = [0,90;45,64.7605981793211;135,64.7605981793211;-135,64.7605981793211;-45,64.7605981793211;0,45;90,45;180,45;-90,45;45,35.2643896827547;135,35.2643896827547;-135,35.2643896827547;-45,35.2643896827547;18.4349488229220,17.5484006137923;71.5650511770780,17.5484006137923;108.434948822922,17.5484006137923;161.565051177078,17.5484006137923;-161.565051177078,17.5484006137923;-108.434948822922,17.5484006137923;-71.5650511770780,17.5484006137923;-18.4349488229220,17.5484006137923;0,0;45,0;90,0;135,0;180,0;-135,0;-90,0;-45,0;18.4349488229220,-17.5484006137923;71.5650511770780,-17.5484006137923;108.434948822922,-17.5484006137923;161.565051177078,-17.5484006137923;-161.565051177078,-17.5484006137923;-108.434948822922,-17.5484006137923;-71.5650511770780,-17.5484006137923;-18.4349488229220,-17.5484006137923;45,-35.2643896827547;135,-35.2643896827547;-135,-35.2643896827547;-45,-35.2643896827547;0,-45;90,-45;180,-45;-90,-45;45,-64.7605981793211;135,-64.7605981793211;-135,-64.7605981793211;-45,-64.7605981793211;0,-90];
        az = azel_rig(lspk,1);
        el = azel_rig(lspk,2);
    end
end

function irBank = extractXrIRs(subjectdir,inv_pol)
    sweepdir = [subjectdir 'sweeps/'];
    [y_inv_sweep, Fs] = audioread([sweepdir 'ZZ_inv_sweep.wav']);

    filelist = dir([sweepdir '*_azi_*.wav']);
    filelist = [filelist; dir([sweepdir '00_reference.wav'])];
    
    sweepBank = struct;
    for i = 1:length(filelist)
        [sweepBank(i).y, sweepBank(i).Fs] = audioread([filelist(i).folder '/' filelist(i).name]);
        sweepBank(i).name = filelist(i).name;

        if inv_pol
            sweepBank(i).y = -sweepBank(i).y;
        end
    end

    irBank = struct;
    for i = 1:length(sweepBank)
        irBank(i).azimuth = str2double(extractBefore(extractAfter(sweepBank(i).name,'_azi_'),'_ele_'));
        irBank(i).elevation = str2double(extractBefore(extractAfter(sweepBank(i).name,'_ele_'),'_dist_'));
        irBank(i).distance = str2double(extractBefore(extractAfter(sweepBank(i).name,'_dist_'),'.wav'));
        irBank(i).Fs = sweepBank(i).Fs;
        irBank(i).name = sweepBank(i).name;
        if contains(sweepBank(i).name,'reference','IgnoreCase',true)
            irBank(i).ref = true;
        else
            irBank(i).ref = false;
        end
        irBank(i).lspk = 1;
        
        fullIr = deconvolve(sweepBank(i).y,y_inv_sweep);
        
        % skip the first half of the ir
        cut_start = floor(size(fullIr,1)/2);
        cut_end = cut_start + 8191;
        irBank(i).fullIR = fullIr(cut_start:cut_end,:);
        
    end
end

function ir = deconvolve(sweep,inv_sweep)
    N1 = size(sweep,1) + size(inv_sweep,1);
    ir(:,1) = ifft((fft(sweep(:,1),N1)).*fft(inv_sweep,N1));
    
    if size(sweep,2) == 2
        ir(:,2) = ifft((fft(sweep(:,2),N1)).*fft(inv_sweep,N1));
    end
end

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

function irBank = adjustGains(irBank,ear_dB,ref_dB)
    ear_gain = 10^(ear_dB/20);
    ref_gain = 10^(ref_dB/20);
    
    for i = 1:length(irBank)
        if irBank(i).ref
            irBank(i).fullIR = irBank(i).fullIR * ref_gain;
        else
            irBank(i).fullIR = irBank(i).fullIR * ear_gain;
        end
    end
end

function plotAndSave(subjectdir,irBank)
    figure('Name','Measured IRs','NumberTitle','off','WindowStyle','docked');
    hold on
    % xlim([0 5000])
    for i = 1:length(irBank)
        plot(irBank(i).fullIR(:,1), '-g')
        if size(irBank(i).fullIR,2) == 2
            plot(irBank(i).fullIR(:,2), '-r')
        end
    end

    save([subjectdir 'irBank.mat'], 'irBank')
end


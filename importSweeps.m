close all
clear

% subjectdir = 'data/20201217-122pt-2.5m-dayton_vt/'; inv_pol = true;
% subjectdir = 'data/20201217-122pt-2.5m-canford_vt/'; inv_pol = true;
% subjectdir = 'data/20211012-q2_tr/'; inv_pol = true;
% subjectdir = 'data/20211105-A-Jan/'; inv_pol = true;
subjectdir = 'data/20211126-TR/'; inv_pol = false;
% subjectdir = 'data/20211126-Gavin/'; inv_pol = false;

sweepdir = [subjectdir 'sweeps/'];
[y_inv_sweep, Fs] = audioread([sweepdir 'ZZ_inv_sweep.wav']);

filelist = dir([sweepdir '*_azi_*.wav']);
filelist = [filelist; dir([sweepdir '00_reference.wav'])];

for i = 1:length(filelist)
    [sweepBank(i).y, sweepBank(i).Fs] = audioread([filelist(i).folder '/' filelist(i).name]);
    sweepBank(i).name = filelist(i).name;
    
    if inv_pol
        sweepBank(i).y = -sweepBank(i).y;
    end
end

for i = 1:length(sweepBank)
    irBank(i).azimuth = str2double(extractBefore(extractAfter(sweepBank(i).name,'_azi_'),'_ele_'));
    irBank(i).elevation = str2double(extractBefore(extractAfter(sweepBank(i).name,'_ele_'),'_dist_'));
    irBank(i).distance = str2double(extractBefore(extractAfter(sweepBank(i).name,'_dist_'),'.wav'));
    irBank(i).Fs = sweepBank(i).Fs;
    irBank(i).name = sweepBank(i).name;
    
    N1 = size(sweepBank(i).y,1) + size(y_inv_sweep,1);
    fullIr = ifft((fft(sweepBank(i).y(:,1),N1)).*fft(y_inv_sweep,N1));
    irLength = size(fullIr,1);
    
    % skip the first half of the ir and cut some silence
%     cut_start = floor(irLength/2) + 2800;
%     cut_end = cut_start + 4095;
    
    cut_start = floor(irLength/2);
    cut_end = cut_start + 8191;
    
    irBank(i).fullIrLeft = fullIr(cut_start:cut_end);
    
    if size(sweepBank(i).y,2) == 2
        fullIr = ifft((fft(sweepBank(i).y(:,2),N1)).*fft(y_inv_sweep,N1));
        irBank(i).fullIrRight = fullIr(cut_start:cut_end);
    end
    


end

%% plot IRs
figure('Name','Measured IRs','NumberTitle','off','WindowStyle','docked');
hold on
% xlim([0 5000])
for i = 1:length(irBank)
    plot(irBank(i).fullIrLeft, '-g')
    plot(irBank(i).fullIrRight, '-r')
end

save([subjectdir 'irBank.mat'], 'irBank')


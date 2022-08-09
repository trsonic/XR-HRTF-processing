close all
clear

%% add paths and initialize SOFA API
addpath('../API_MO/SOFAtoolbox/')
SOFAstart;

sofafile = 'data/20220807-XR-SUBJ001B/xr-hrtf-interp-raw.sofa';

%% load sofa file
hrtf = SOFAload(sofafile);
hrirBank = [];
apparentSourceVector = SOFAcalculateAPV(hrtf); % Calculate the source position from a listener point of view
for i = 1:length(apparentSourceVector)
    hrirBank(i).azimuth = apparentSourceVector(i,1);
    hrirBank(i).elevation = apparentSourceVector(i,2);
    hrirBank(i).left_hrir = squeeze(hrtf.Data.IR(i,1,:))';
    hrirBank(i).right_hrir = squeeze(hrtf.Data.IR(i,2,:))';
    hrirBank(i).Fs = hrtf.Data.SamplingRate;        
end

% remove duplicated directions (just in case)
[~, idx] = unique([[hrirBank.azimuth]',[hrirBank.elevation]'],'rows');
hrirBank = hrirBank(idx);

%% plot surface plots
plotHMFmags(hrirBank)




close all
clear

%% add paths and initialize SOFA API
addpath('../API_MO/API_MO/')
SOFAstart;

sofafile = 'data/20211126-XR-TR/xr-hrtf-interp-raw.sofa';
% sofafile = 'data/20220223-XR-TR_median/xr-hrtf-raw.sofa';


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




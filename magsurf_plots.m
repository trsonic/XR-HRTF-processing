close all
clear

%% add paths and initialize SOFA API
addpath('../API_MO/API_MO/')
SOFAstart;

sofafile = 'data/20201217-122pt-2.5m-dayton_vt/xr_122pt.sofa';

%% load sofa file
hrtf = SOFAload(sofafile);
hrirBank = [];
apparentSourceVector = SOFAcalculateAPV(hrtf); % Calculate the source position from a listener point of view
for i = 1:length(apparentSourceVector)
    hrirBank(i).azimuth = apparentSourceVector(i,1);
    hrirBank(i).elevation = apparentSourceVector(i,2);
%     hrirBank(i).rawhrir = squeeze(hrtf.Data.IR(i,:,:))';
    hrirBank(i).left_hrir = squeeze(hrtf.Data.IR(i,1,:))';
    hrirBank(i).right_hrir = squeeze(hrtf.Data.IR(i,2,:))';
    hrirBank(i).Fs = hrtf.Data.SamplingRate;        
end

% remove duplicated directions (just in case)
[~, idx] = unique([[hrirBank.azimuth]',[hrirBank.elevation]'],'rows');
hrirBank = hrirBank(idx);

%% plot surface plots
figure('Name','Surf plots','NumberTitle','off','WindowStyle','docked')
tiledlayout(3,2)

% horizontal plane
az = -180:1:180;
el = 0;
left_hrirs = [];
right_hrirs = [];
for i = 1:length(az)
    dist = distance([hrirBank.elevation], [hrirBank.azimuth], el, az(i));
    [~, idx] = min(dist);
    left_hrirs(:,i) = hrirBank(idx).left_hrir;
    right_hrirs(:,i) = hrirBank(idx).right_hrir;
end
nexttile
magsurf(hrirBank(idx).Fs,az,left_hrirs,'Azimuth (deg)','Horizontal plane - left ear',[-60 20])
nexttile
magsurf(hrirBank(idx).Fs,az,right_hrirs,'Azimuth (deg)','Horizontal plane - right ear',[-60 20])

% median plane
az = 0;
el = -180:1:180;
for i = 1:length(el)
    dist = distance([hrirBank.elevation], [hrirBank.azimuth], el(i), az);
    [~, idx] = min(dist);
    left_hrirs(:,i) = hrirBank(idx).left_hrir;
    right_hrirs(:,i) = hrirBank(idx).right_hrir;
end
nexttile
magsurf(hrirBank(idx).Fs,el,left_hrirs,'Elevation (deg)','Median plane - left ear',[-60 20])
nexttile
magsurf(hrirBank(idx).Fs,el,right_hrirs,'Elevation (deg)','Median plane - right ear',[-60 20])

% frontal plane
az = 90;
el = -180:1:180;
for i = 1:length(el)
    dist = distance([hrirBank.elevation], [hrirBank.azimuth], el(i), az);
    [~, idx] = min(dist);
    left_hrirs(:,i) = hrirBank(idx).left_hrir;
    right_hrirs(:,i) = hrirBank(idx).right_hrir;
end
nexttile
magsurf(hrirBank(idx).Fs,el,left_hrirs,'Elevation (deg)','Frontal plane - left ear',[-60 20])
nexttile
magsurf(hrirBank(idx).Fs,el,right_hrirs,'Elevation (deg)','Frontal plane - right ear',[-60 20])

function magsurf(Fs,ang,hrir,ytit,tit,zlimit)
    Nfft = 1024;
    f = (Fs/Nfft:Fs/Nfft:Fs)';
    magnitude = abs(fft(hrir,Nfft));
    magnitude = 20*log10(magnitude);
    s = pcolor(f,ang,magnitude');
    s.EdgeColor = 'none';
    hold on
    yline(-90,'--')
    yline(0,'--')
    yline(90,'--')
    set(gca,'xscale','log')
    xlim([500 24000])
    xticks([1000 2000 4000 8000 16000])
    ylim([-180 180])
    zlim(zlimit)
    caxis(zlimit)
    xlabel('Frequency (Hz)')
    ylabel(ytit)
    title(tit)
    colorbar
end

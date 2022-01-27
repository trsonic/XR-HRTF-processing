close all
clear

addpath('tools/')
% addpath('../XR-HRTF-processing/tools/VoronoiSphere/')

%% create a huge IRbank with all windowed KEMAR measurements
dirlist = dir('data/*-RIG-*');
IRbank = [];
for i = 1:length(dirlist)
    subjectdir = [dirlist(i).folder '/' dirlist(i).name '/'];
    load([subjectdir 'irBankProcessed.mat'])

    for j = 1:size(irBank,2)
        for ch = 1:size(irBank(i).winIR,2)
           IRbank(end+1).name = irBank(j).name;
           IRbank(end).type = extractBefore(extractAfter(IRbank(end).name,'-'),'-');
           IRbank(end).hmd = extractBefore(extractAfter(IRbank(end).name,'deg-'),'-');
           IRbank(end).deg = extractAfter(extractBefore(IRbank(end).name,'deg'),'KEMAR-');
           IRbank(end).lspk = irBank(j).lspk;
           IRbank(end).ch = ch;
           if ch == 1
               IRbank(end).azimuth = irBank(j).azimuth;
               IRbank(end).ITD = irBank(j).ITD;
               IRbank(end).dly = irBank(j).dlyL;
           else
               IRbank(end).azimuth = -irBank(j).azimuth;
               IRbank(end).ITD = -irBank(j).ITD;
               IRbank(end).dly = irBank(j).dlyR;
           end
           IRbank(end).elevation = irBank(j).elevation;
           IRbank(end).winIR = irBank(j).winIR(:,ch);
           IRbank(end).Fs = irBank(j).Fs;
        end
    end
end

%% get magnitude vectors
for i = 1:size(IRbank,2)
    [IRbank(i).f,IRbank(i).mag] = getMagnitude(IRbank(i).winIR,IRbank(i).Fs,'lin');
end


%% create single IR banks
IRb_NOHMD = IRbank(getIdx(IRbank,'hmd','NOHMD'));
IRb_NOHMD2 = IRbank(getIdx(IRbank,'hmd','NOHMD2'));
IRb_Q2HMD = IRbank(getIdx(IRbank,'hmd','Q2HMD'));
IRb_MLHMD = IRbank(getIdx(IRbank,'hmd','MLHMD'));
IRb_TRHAT = IRbank(getIdx(IRbank,'hmd','TRHAT'));

%% calculate TOD difference for Quest
model=struct;
for i = 1:size(IRb_Q2HMD,2)
    model(i).azi = IRb_Q2HMD(i).azimuth;
    model(i).ele = IRb_Q2HMD(i).elevation;
    model(i).dtoa_diff = (IRb_Q2HMD(i).dly-(IRb_NOHMD(i).dly+IRb_NOHMD2(i).dly)/2)';
%     model(i).mag_diff = (IRb_Q2HMD(i).mag-(IRb_NOHMD(i).mag+IRb_NOHMD2(i).mag)/2)'; %%%%%
    model(i).mag_diff = (IRb_Q2HMD(i).mag ./ sqrt(IRb_NOHMD(i).mag.^2+IRb_NOHMD2(i).mag.^2))';    
    model(i).f = IRb_Q2HMD(i).f;
end


ls = getLebedevSphere(4334);
dirs = [];
[dirs(:,1), dirs(:,2), ~] = cart2sph(ls.x,ls.y,ls.z);
dirs = rad2deg(dirs);
% [model.azi] = dirs(:,1);
% model.ele = dirs(:,2)';

% gaussian window
sigma = 5;
pd = makedist('Normal','mu',0,'sigma',sigma);
% plot(pdf(pd,linspace(0,180,180)))


model_interp = struct;
for i = 1:size(dirs,1)
    model_interp(i).az = dirs(i,1);
    model_interp(i).el = dirs(i,2);
    
    dist = distance([model.ele],[model.azi],dirs(i,2),dirs(i,1));
    [val,idx] = min(dist);
    
    weights = pdf(pd,dist);
    weights = weights / sum(weights);
    model_interp(i).dtoa_diff = sum([model.dtoa_diff] .* weights);
    
    % calculate magnitude difference
    mag_diff = zeros(1,size(model(1).mag_diff,2));

    for j = 1:size(model,2)
        mag_diff = mag_diff + (model(j).mag_diff.^2) * weights(j);
    end
    f = model(1).f;
    mag_diff = sqrt(mag_diff);
    model_interp(i).f = f;
    model_interp(i).mag_diff = mag_diff;
    mag_diff_mf = 20*log10(rms(mag_diff(f >= 2000 & f <= 8000)));
    model_interp(i).mag_diff_mf = mag_diff_mf;
end




% %% plot quest TOA difference interpolated
% figure('Name','TOA difference','NumberTitle','off','WindowStyle','docked');
% % tiledlayout(2,2)
% lim = [-20 80];
% nexttile
% hold on
% title('quest TOA difference interpolated')
% plotAzEl([model_interp.az],[model_interp.el],[model_interp.dtoa_diff],lim)
% 
% 
% %% plot quest magnitude difference interpolated
% figure('Name','mag difference','NumberTitle','off','WindowStyle','docked');
% % tiledlayout(2,2)
% lim = [-12 12];
% nexttile
% hold on
% title('quest magnitude difference interpolated')
% plotAzEl([model_interp.az],[model_interp.el],[model_interp.mag_diff_mf],lim)


%% plot quest magnitude difference interpolated
figure('Name','mag difference','NumberTitle','off','WindowStyle','docked');
hold on

for i = 1:size(model_interp,2)
	plot(model_interp(i).f,20*log10(model_interp(i).mag_diff));
end

set(gca,'xscale','log')
grid on
xlim([10 24000]);



%% plot TOA difference
figure('Name','TOA difference','NumberTitle','off','WindowStyle','docked');
tiledlayout(2,2)
lim = [-20 80];
nexttile
hold on
title('repeated')
diff = [IRb_NOHMD2.dly]-[IRb_NOHMD.dly];
plotAzEl([IRb_NOHMD.azimuth],[IRb_NOHMD.elevation],diff,lim)

nexttile
hold on
title('Quest 2')
diff = [IRb_Q2HMD.dly]-([IRb_NOHMD.dly]+[IRb_NOHMD2.dly])/2;
disp(['maximum TOA difference: ' num2str(max(abs(diff))) ' us'])
plotAzEl([IRb_NOHMD.azimuth],[IRb_NOHMD.elevation],diff,lim)

nexttile
hold on
title('Magic Leap One')
diff = [IRb_MLHMD.dly]-([IRb_NOHMD.dly]+[IRb_NOHMD2.dly])/2;
plotAzEl([IRb_NOHMD.azimuth],[IRb_NOHMD.elevation],diff,lim)

nexttile
hold on
title('Woolen hat')
diff = [IRb_TRHAT.dly]-([IRb_NOHMD.dly]+[IRb_NOHMD2.dly])/2;
plotAzEl([IRb_NOHMD.azimuth],[IRb_NOHMD.elevation],diff,lim)

% %% plot ITD
% figure('Name','ITD','NumberTitle','off','WindowStyle','docked');
% tiledlayout(2,2)
% lim = [-800 800];
% nexttile
% hold on
% plotAzEl([IRb_NOHMD.azimuth],[IRb_NOHMD.elevation],[IRb_NOHMD.ITD], lim)
% title('No HMD')
% nexttile
% hold on
% plotAzEl([IRb_NOHMD.azimuth],[IRb_NOHMD.elevation],[IRb_Q2HMD.ITD], lim)
% title('Quest 2')
% nexttile
% hold on
% plotAzEl([IRb_NOHMD.azimuth],[IRb_NOHMD.elevation],[IRb_MLHMD.ITD], lim)
% title('Magic Leap One')
% nexttile
% hold on
% plotAzEl([IRb_NOHMD.azimuth],[IRb_NOHMD.elevation],[IRb_TRHAT.ITD], lim)
% title('Woolen hat')

%% plot ITD difference
figure('Name','ITD difference','NumberTitle','off','WindowStyle','docked');
tiledlayout(2,2)
lim = [-100 100];
nexttile
hold on
plotAzEl([IRb_NOHMD.azimuth],[IRb_NOHMD.elevation],[IRb_NOHMD2.ITD]-[IRb_NOHMD.ITD], lim)
title('repeated')
nexttile
hold on
title('Quest 2')
diff = [IRb_Q2HMD.ITD]-([IRb_NOHMD.ITD]+[IRb_NOHMD2.ITD])/2;
disp(['maximum ITD difference: ' num2str(max(abs(diff))) ' us'])
plotAzEl([IRb_NOHMD.azimuth],[IRb_NOHMD.elevation],diff,lim)

nexttile
hold on
plotAzEl([IRb_NOHMD.azimuth],[IRb_NOHMD.elevation],[IRb_MLHMD.ITD]-[IRb_NOHMD.ITD], lim)
title('Magic Leap One')
nexttile
hold on
plotAzEl([IRb_NOHMD.azimuth],[IRb_NOHMD.elevation],[IRb_TRHAT.ITD]-[IRb_NOHMD.ITD], lim)
title('Woolen hat')

%% plot ITD difference (above JND)
figure('Name','ITD difference','NumberTitle','off','WindowStyle','docked');
tiledlayout(2,2)
lim = [-100 100];
nexttile
hold on
val = clipValues([IRb_NOHMD2.ITD]-[IRb_NOHMD.ITD]);
plotAzEl([IRb_NOHMD.azimuth],[IRb_NOHMD.elevation],val, lim)
title('repeated')
nexttile
hold on
val = mean([[IRb_Q2HMD.ITD]-[IRb_NOHMD.ITD]; [IRb_Q2HMD.ITD]-[IRb_NOHMD2.ITD]],1);
val = clipValues(val);
plotAzEl([IRb_NOHMD.azimuth],[IRb_NOHMD.elevation],val, lim)
title('Quest 2')
nexttile
hold on
val = clipValues([IRb_TRHAT.ITD]-[IRb_NOHMD.ITD]);
val = mean([[IRb_MLHMD.ITD]-[IRb_NOHMD.ITD]; [IRb_MLHMD.ITD]-[IRb_NOHMD2.ITD]],1);
val = clipValues(val);
plotAzEl([IRb_NOHMD.azimuth],[IRb_NOHMD.elevation],val, lim)
title('Magic Leap One')
nexttile
hold on
val = mean([[IRb_TRHAT.ITD]-[IRb_NOHMD.ITD]; [IRb_TRHAT.ITD]-[IRb_NOHMD2.ITD]],1);
val = clipValues(val);
plotAzEl([IRb_NOHMD.azimuth],[IRb_NOHMD.elevation],val, lim)
title('Woolen hat')


%% functions
function val = clipValues(val)
    val(abs(val) < 15) = 0;
end

function idx = getIdx(IRbank,cat,key)
    idx = [];
    for i = 1:length(IRbank)
        if strcmp('hmd', cat)
            if strcmp(IRbank(i).hmd, key)
                idx = [idx; i];
            end
        end
    end
end

function plotAzEl(az,el,val,lim)
    azimuth = -180:1:180;
    elevation = -90:1:90;

    for i = 1:length(azimuth)
        for j = 1:length(elevation)
            dist = distance(elevation(j),azimuth(i),el,az);
            idx = find(dist == min(dist));
            value(j,i) = mean(val(idx));
%             if length(idx) > 2
%                 disp(val(idx));
%             end
        end    
    end

    s = pcolor(azimuth,elevation,value);
    s.EdgeColor = 'none';
    xlim([-180 180])
    ylim([-90 90])
    zlim(lim)
    caxis(lim)
    xlabel('Azimuth (deg)')
    ylabel('Elevation (deg)')
    colorbar
end

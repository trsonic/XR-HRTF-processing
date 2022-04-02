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

%% ERB stuff
erb = 1:40;
Q = 9.265;
L = 24.7;
gamma = 1;
for i = erb
    fc(i) = Q * L * (exp(i*gamma/Q) - 1);
    fbw(i) = gamma * L * exp(i*gamma/Q);
end

% figure
% hold on
% plot(erb,fc)
% plot(erb,erb2hz(erb))


%% get magnitude vectors and ERB
for i = 1:size(IRbank,2)
    [f,mag] = getMagnitude(IRbank(i).winIR,IRbank(i).Fs,'lin');
    
    IRbank(i).f = f;
    IRbank(i).mag = mag;
    
    for j = erb
       erb_pwr(j) = sum(mag(f >= fc(j)-fbw(j)/2 & f < fc(j)+fbw(j)/2).^2);
    end
    IRbank(i).erb_pwr = erb_pwr;
end


%% create single IR banks
IRb_NOHMD = IRbank(getIdx(IRbank,'hmd','NOHMD'));
IRb_NOHMD2 = IRbank(getIdx(IRbank,'hmd','NOHMD2'));
IRb_Q2HMD = IRbank(getIdx(IRbank,'hmd','Q2HMD'));
IRb_MLHMD = IRbank(getIdx(IRbank,'hmd','MLHMD'));
IRb_TRHAT = IRbank(getIdx(IRbank,'hmd','TRHAT'));

%% calculate TOA and magnitude difference for Quest
model=struct;
for i = 1:size(IRb_Q2HMD,2)
    model(i).azi = IRb_Q2HMD(i).azimuth;
    model(i).ele = IRb_Q2HMD(i).elevation;
    model(i).dtoa_diff = (IRb_Q2HMD(i).dly-(IRb_NOHMD(i).dly+IRb_NOHMD2(i).dly)/2)';
%     model(i).mag_diff = (IRb_Q2HMD(i).mag-(IRb_NOHMD(i).mag+IRb_NOHMD2(i).mag)/2)';
%     model(i).mag_diff = (IRb_Q2HMD(i).mag ./ sqrt(IRb_NOHMD(i).mag.^2+IRb_NOHMD2(i).mag.^2))';
    model(i).mag_diff = IRb_Q2HMD(i).mag ./ sqrt(0.5*IRb_NOHMD(i).mag.^2+0.5*IRb_NOHMD2(i).mag.^2);    
    
    for j = erb
        model(i).sd(j) =  10 * log10(IRb_Q2HMD(i).erb_pwr(j) ./ sqrt(0.5*IRb_NOHMD(i).erb_pwr(j)^2+0.5*IRb_NOHMD2(i).erb_pwr(j)^2));
    end
    model(i).f = IRb_Q2HMD(i).f;
end

%% plot analyzed directions
dirs = unique([[model.azi]' [model.ele]'],'rows');
idx_to_remove = [];
for i = 1:size(dirs,1)
    dist = distance(dirs(:,2),dirs(:,1),dirs(i,2),dirs(i,1));
    idx = find(dist < 0.1);
    
    if length(idx) > 1
       disp(idx)
       idx_to_remove = [idx_to_remove; idx(2:end)];
    end
    
end

idx_to_remove = unique(idx_to_remove);
dirs2 = dirs;
dirs2(idx_to_remove,:)=[];

figure('Name','directions mollweide','NumberTitle','off','WindowStyle','docked');
hold on


axesm('MapProjection','mollweid','MapLatLimit',[-90 90],'Gcolor','black','GLineWidth',1.0,'MLineLocation',45,'PLineLocation',30); 
axis off
gridm('on')
% scatterm(dirs(:,2),dirs(:,1))
scatterm(dirs2(:,2),dirs2(:,1),120,'k')

% labels
fsize = 12;
fcolor = 'black';
vertshift = -5;
horishift = 5;
textm(vertshift,horishift-135,'-135','color',fcolor,'fontsize',fsize);
textm(vertshift,horishift-90,'-90','color',fcolor,'fontsize',fsize);
textm(vertshift,horishift-45,'-45','color',fcolor,'fontsize',fsize);
textm(vertshift,horishift+0,'0','color',fcolor,'fontsize',fsize);
textm(vertshift,horishift+45,'45','color',fcolor,'fontsize',fsize);
textm(vertshift,horishift+90,'90','color',fcolor,'fontsize',fsize);
textm(vertshift,horishift+135,'135','color',fcolor,'fontsize',fsize);
textm(-vertshift-60,horishift+0,'-60','color',fcolor,'fontsize',fsize);
textm(-vertshift-30,horishift+0,'-30','color',fcolor,'fontsize',fsize);
textm(vertshift+30,horishift+0,'30','color',fcolor,'fontsize',fsize);
textm(vertshift+60,horishift+0,'60','color',fcolor,'fontsize',fsize);

% save figure
figlen = 8;
width = 4*figlen;
height = 3*figlen;
set(gcf,'Units','centimeters','PaperPosition',[0 0 width height],'PaperSize',[width height]);
saveas(gcf,'data/hmdpert_output/analyzed_dirs.png')

%% interpolate these differences using gaussian smoothing
% ls = getLebedevSphere(350);
ls = getLebedevSphere(4334);
% ls = getLebedevSphere(50);
dirs = [];
[dirs(:,1), dirs(:,2), ~] = cart2sph(ls.x,ls.y,ls.z);
dirs = rad2deg(dirs);

% gaussian window
sigma = 5;
pd = makedist('Normal','mu',0,'sigma',sigma);
figure('Name','smoothing window','NumberTitle','off','WindowStyle','docked');

distvec = 0:0.25:180;
plot(distvec,pdf(pd,distvec))
% ylim([0 1])
xlim([0 60])
ylabel('Weighting')
xlabel('Angular Distance (deg)')

% create struct with interpolated values
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
    sd = zeros(1,size(model(1).sd,2));

    for j = 1:size(model,2)
        mag_diff = mag_diff + (model(j).mag_diff.^2) * weights(j);
        sd = sd + model(j).sd * weights(j);
    end
    f = model(1).f;
    mag_diff = sqrt(mag_diff);
    
    model_interp(i).f = f;
    model_interp(i).mag_diff = mag_diff;
    model_interp(i).mag_diff_lf = rms(10*log10(mag_diff(f >= 100 & f < 1000))); % low
    model_interp(i).mag_diff_mf = rms(10*log10(mag_diff(f >= 1000 & f < 5000))); % mid
    model_interp(i).mag_diff_hf = rms(10*log10(mag_diff(f >= 5000 & f < 16000))); % high
    
%     model_interp(i).sd_lf = sqrt(mean(sd(fc >= 100 & fc < 1000).^2));       % low
%     model_interp(i).sd_mf = sqrt(mean(sd(fc >= 1000 & fc < 5000).^2));      % mid
%     model_interp(i).sd_hf = sqrt(mean(sd(fc >= 5000 & fc < 16000).^2));     % high

    model_interp(i).sd_lf = mean(sd(fc >= 100 & fc < 1000),'omitnan');       % low
    model_interp(i).sd_mf = mean(sd(fc >= 1000 & fc < 5000),'omitnan');      % mid
    model_interp(i).sd_hf = mean(sd(fc >= 5000 & fc < 16000),'omitnan');     % high


end

%% plot quest TOA difference interpolated
figure('Name','quest TOA difference interpolated','NumberTitle','off','WindowStyle','docked');
% tiledlayout(1,2)
lim = [-20 60];
% nexttile
hold on
% title('quest TOA difference interpolated')
plotAzElM([model_interp.az],[model_interp.el],[model_interp.dtoa_diff],lim,'Time-of-arrival difference (us)')

% nexttile
% hold on
% title('quest TOA difference interpolated')
% plotAzEl([model_interp.az],[model_interp.el],[model_interp.dtoa_diff],lim)

% save figure
figlen = 8;
width = 4*figlen;
height = 3*figlen;
set(gcf,'Units','centimeters','PaperPosition',[0 0 width height],'PaperSize',[width height]);
saveas(gcf,'data/hmdpert_output/toa_difference.png')


%% plot quest magnitude difference interpolated
figure('Name','quest magnitude difference interpolated','NumberTitle','off','WindowStyle','docked');
t = tiledlayout(3,1);
t.TileSpacing = 'compact';
t.Padding = 'compact';
nexttile
lim = [-8 8];
hold on
plotAzElM([model_interp.az],[model_interp.el],[model_interp.sd_lf],lim,'Magnitude difference (dB)')

nexttile
lim = [-8 8];
hold on
plotAzElM([model_interp.az],[model_interp.el],[model_interp.sd_mf],lim,'Magnitude difference (dB)')

nexttile
lim = [-8 8];
hold on
plotAzElM([model_interp.az],[model_interp.el],[model_interp.sd_hf],lim,'Magnitude difference (dB)')

% save figure
figlen = 8;
width = 4*figlen;
height = 7*figlen;
set(gcf,'Units','centimeters','PaperPosition',[0 0 width height],'PaperSize',[width height]);
saveas(gcf,'data/hmdpert_output/mag_difference.png')  


%% plot quest magnitude difference interpolated
figure('Name','quest magnitude difference interpolated','NumberTitle','off','WindowStyle','docked');
t = tiledlayout(3,1);
t.TileSpacing = 'compact';
t.Padding = 'compact';
nexttile
lim = [-12 12];
hold on
plotAzElM([model_interp.az],[model_interp.el],[model_interp.mag_diff_lf],lim,'Magnitude difference (dB)')

nexttile
lim = [-12 12];
hold on
plotAzElM([model_interp.az],[model_interp.el],[model_interp.mag_diff_mf],lim,'Magnitude difference (dB)')

nexttile
lim = [-12 12];
hold on
plotAzElM([model_interp.az],[model_interp.el],[model_interp.mag_diff_hf],lim,'Magnitude difference (dB)')

% save figure
figlen = 8;
width = 4*figlen;
height = 7*figlen;
set(gcf,'Units','centimeters','PaperPosition',[0 0 width height],'PaperSize',[width height]);
saveas(gcf,'data/hmdpert_output/mag_difference.png')  

%% CREATE INVERSE FILTERS

%% plot quest magnitude difference interpolated
figure('Name','mag difference','NumberTitle','off','WindowStyle','docked');

for i = 1:size(model_interp,2)
    f = model_interp(i).f;
    H = model_interp(i).mag_diff;
    
    H_sm = smoothSpectrum(H(1:end/2),f(1:end/2),3);
    H_sm = [H_sm; fliplr(H_sm')'];
    
    iH = conj(H_sm)./(conj(H_sm).*H_sm);
    invh = circshift(ifft(iH,'symmetric'),length(iH)/2);
    invh = minph(invh);
    
    subplot(2,1,1)
    hold on

% 	plot(f,20*log10(H))
%     plot(f,20*log10(H_sm))
    plot(f,20*log10(iH))
    
    set(gca,'xscale','log')
    grid on
    xlim([10 24000])
    ylim([-12 12])
    
    subplot(2,1,2)
    hold on
    plot(invh)
    
    model_interp(i).invh = invh;
%     subplot(3,1,3)
%     freqz(invh,0:100:24000,48000)
    
end

save('data/hmdpert_output/model_interp.mat', 'model_interp')

%% RAW MEASUREMENT DATA (NOT USED)
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
%     val(abs(val) < 15) = 0;
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

function plotAzElM(az,el,val,lim,clbl)
    azimuth = -180:1:180;
    elevation = -90:1:90;

    for i = 1:length(azimuth)
        for j = 1:length(elevation)
            dist = distance(elevation(j),azimuth(i),el,az);
            idx = find(dist == min(dist));
            value(j,i) = mean(val(idx));
        end    
    end
    
    axesm('MapProjection','mollweid','MapLatLimit',[-90 90],'Gcolor','black','GLineWidth',1.0,'MLineLocation',45,'PLineLocation',30); 
    axis off
    gridm('on')

    % plot data
    s = pcolorm(elevation,azimuth,value);
    s.EdgeColor = 'none';

    % color bar
    caxis(lim)
%     c=colorbar('location','southoutside','position',[0.2 0.0 0.6 0.05],'box','on','color',[0 0 0]); %[left, bottom, width, height]
    c=colorbar('location','southoutside');
    c.Label.String=clbl;
    c.Label.FontSize=14;
    c.FontSize=14;
%     c.Limits=[0 500];
%     c.Ticks=0:50:500;

    % labels
    fsize = 14;
    fcolor = 'black';
    vertshift = -5;
    horishift = 5;
    textm(vertshift,horishift-135,'-135','color',fcolor,'fontsize',fsize);
    textm(vertshift,horishift-90,'-90','color',fcolor,'fontsize',fsize);
    textm(vertshift,horishift-45,'-45','color',fcolor,'fontsize',fsize);
    textm(vertshift,horishift+0,'0','color',fcolor,'fontsize',fsize);
    textm(vertshift,horishift+45,'45','color',fcolor,'fontsize',fsize);
    textm(vertshift,horishift+90,'90','color',fcolor,'fontsize',fsize);
    textm(vertshift,horishift+135,'135','color',fcolor,'fontsize',fsize);
    textm(vertshift-60,horishift+0,'-60','color',fcolor,'fontsize',fsize);
    textm(vertshift-30,horishift+0,'-30','color',fcolor,'fontsize',fsize);
    textm(vertshift+30,horishift+0,'30','color',fcolor,'fontsize',fsize);
    textm(vertshift+60,horishift+0,'60','color',fcolor,'fontsize',fsize);
    
    text(-2,1.25,'contralateral','color',fcolor,'fontsize',fsize,'rotation',0,'horizontalalignment','center','verticalalignment','middle');
    text(2,1.25,'ipsilateral','color',fcolor,'fontsize',fsize,'rotation',0,'horizontalalignment','center','verticalalignment','middle');
    
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

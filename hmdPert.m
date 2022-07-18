close all
clear

addpath('tools/')
addpath('tools/VoronoiSphere/')

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
               IRbank(end).ILD = irBank(j).ILD;
           else
               IRbank(end).azimuth = -irBank(j).azimuth;
               IRbank(end).ITD = -irBank(j).ITD;
               IRbank(end).dly = irBank(j).dlyR;
               IRbank(end).ILD = -irBank(j).ILD;
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
    model(i).itd_diff = (IRb_Q2HMD(i).ITD-(IRb_NOHMD(i).ITD+IRb_NOHMD2(i).ITD)/2)';
    model(i).dtoa_diff = (IRb_Q2HMD(i).dly-(IRb_NOHMD(i).dly+IRb_NOHMD2(i).dly)/2)';
    model(i).ild_diff = (IRb_Q2HMD(i).ILD-(IRb_NOHMD(i).ILD+IRb_NOHMD2(i).ILD)/2)';
    model(i).mag_diff = IRb_Q2HMD(i).mag ./ sqrt(0.5*IRb_NOHMD(i).mag.^2+0.5*IRb_NOHMD2(i).mag.^2);    
    
    for j = erb
        model(i).sd(j) =  10 * log10(IRb_Q2HMD(i).erb_pwr(j) ./ sqrt(0.5*IRb_NOHMD(i).erb_pwr(j)^2+0.5*IRb_NOHMD2(i).erb_pwr(j)^2));
    end
    model(i).f = IRb_Q2HMD(i).f;
end

%% calculate spectral difference
for i = 1:size(model,2)
    model(i).sd_lf = rms(model(i).sd(fc >= 100 & fc < 1000),'omitnan');       % low
    model(i).sd_mf = rms(model(i).sd(fc >= 1000 & fc < 5000),'omitnan');      % mid
    model(i).sd_hf = rms(model(i).sd(fc >= 5000 & fc < 16000),'omitnan');     % high
end



%% plot quest ITD Error
figure('Name','quest ITD Error','NumberTitle','off','WindowStyle','docked');
% tiledlayout(1,2)
lim = [0 80];
% nexttile
hold on
plotAzElM([model.azi],[model.ele],abs([model.itd_diff]),lim,'ITD Error (μs)','','')

% save figure
figlen = 8;
width = 4*figlen;
height = 3*figlen;
set(gcf,'Units','centimeters','PaperPosition',[0 0 width height],'PaperSize',[width height]);
saveas(gcf,'data/hmdpert_output/itd_error.png')

%% plot quest ILD Error
figure('Name','quest ILD Error','NumberTitle','off','WindowStyle','docked');
% tiledlayout(1,2)
lim = [0 4];
% nexttile
hold on
plotAzElM([model.azi],[model.ele],abs([model.ild_diff]),lim,'ILD Error (dB)','','')

% [maxv,idx] = max(abs([model.ild_diff]));


% save figure
figlen = 8;
width = 4*figlen;
height = 3*figlen;
set(gcf,'Units','centimeters','PaperPosition',[0 0 width height],'PaperSize',[width height]);
saveas(gcf,'data/hmdpert_output/ild_error.png')

% %% plot quest TOA difference
% figure('Name','quest TOA difference','NumberTitle','off','WindowStyle','docked');
% % tiledlayout(1,2)
% lim = [-20 80];
% % nexttile
% hold on
% plotAzElM([model.azi],[model.ele],[model.dtoa_diff],lim,'Time-of-arrival difference (μs)','contralateral','ipsilateral')
% 
% % save figure
% figlen = 8;
% width = 4*figlen;
% height = 3*figlen;
% set(gcf,'Units','centimeters','PaperPosition',[0 0 width height],'PaperSize',[width height]);
% saveas(gcf,'data/hmdpert_output/toa_difference.png')
% % exportgraphics(gcf,'data/hmdpert_output/toa_difference.png','Resolution',300)


%% plot quest spectral difference
figure('Name','quest spectral difference','NumberTitle','off','WindowStyle','docked');
t = tiledlayout(3,1);
t.TileSpacing = 'compact';
t.Padding = 'compact';
nexttile
lim = [0 8];
hold on
plotAzElM([model.azi],[model.ele],[model.sd_lf],lim,'Spectral difference (dB)','contralateral','ipsilateral')
annotation('textbox', [.055, (3/3)-0.1-0.01, .15, .05], 'string', '0.1 - 1 kHz','fontsize',18,'LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle')

nexttile
lim = [0 8];
hold on
plotAzElM([model.azi],[model.ele],[model.sd_mf],lim,'Spectral difference (dB)','contralateral','ipsilateral')
annotation('textbox', [.055, (2/3)-0.1, .15, .05], 'string', '1 - 5 kHz','fontsize',18,'LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle')

nexttile
lim = [0 8];
hold on
plotAzElM([model.azi],[model.ele],[model.sd_hf],lim,'Spectral difference (dB)','contralateral','ipsilateral')
annotation('textbox', [.055, (1/3)-0.1+0.01, .15, .05], 'string', '5 - 16 kHz','fontsize',18,'LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle')

% save figure
figlen = 8;
width = 4*figlen;
height = 7*figlen;
set(gcf,'Units','centimeters','PaperPosition',[0 0 width height],'PaperSize',[width height]);
saveas(gcf,'data/hmdpert_output/spectral_difference.png')  


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
    
    weights = pdf(pd,dist);
    weights = weights / sum(weights);
    model_interp(i).dtoa_diff = sum([model.dtoa_diff] .* weights);
    model_interp(i).f = model(1).f;
    
    % calculate interpolated magnitude difference
    mag_diff = zeros(1,size(model(1).mag_diff,2));
    for j = 1:size(model,2)
        mag_diff = mag_diff + (model(j).mag_diff.^2) * weights(j);
    end
    model_interp(i).mag_diff = sqrt(mag_diff);

    % calculate interpolated sd
    sd = zeros(1,size(model(1).sd,2));
    for j = 1:size(model,2)
        sd = sd + model(j).sd * weights(j);
    end
    
    model_interp(i).sd_lf = rms(sd(fc >= 100 & fc < 1000),'omitnan');       % low
    model_interp(i).sd_mf = rms(sd(fc >= 1000 & fc < 5000),'omitnan');      % mid
    model_interp(i).sd_hf = rms(sd(fc >= 5000 & fc < 16000),'omitnan');     % high
end

%% plot quest TOA difference interpolated
figure('Name','quest TOA difference interpolated','NumberTitle','off','WindowStyle','docked');
% tiledlayout(1,2)
lim = [-20 80];
% nexttile
hold on
% title('quest TOA difference interpolated')
plotAzElM([model_interp.az],[model_interp.el],[model_interp.dtoa_diff],lim,'Time-of-arrival difference (μs)','contralateral','ipsilateral')

% save figure
figlen = 8;
width = 4*figlen;
height = 3*figlen;
set(gcf,'Units','centimeters','PaperPosition',[0 0 width height],'PaperSize',[width height]);
saveas(gcf,'data/hmdpert_output/toa_difference_interpolated.png')


%% plot quest magnitude difference interpolated
figure('Name','quest magnitude difference interpolated','NumberTitle','off','WindowStyle','docked');
t = tiledlayout(3,1);
t.TileSpacing = 'compact';
t.Padding = 'compact';
nexttile
lim = [-8 8];
hold on
plotAzElM([model_interp.az],[model_interp.el],[model_interp.sd_lf],lim,'Magnitude difference (dB)','contralateral','ipsilateral')
annotation('textbox', [.055, (3/3)-0.1-0.01, .15, .05], 'string', '0.1 - 1 kHz','fontsize',18,'LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle')

nexttile
lim = [-8 8];
hold on
plotAzElM([model_interp.az],[model_interp.el],[model_interp.sd_mf],lim,'Magnitude difference (dB)','contralateral','ipsilateral')
annotation('textbox', [.055, (2/3)-0.1, .15, .05], 'string', '1 - 5 kHz','fontsize',18,'LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle')

nexttile
lim = [-8 8];
hold on
plotAzElM([model_interp.az],[model_interp.el],[model_interp.sd_hf],lim,'Magnitude difference (dB)','contralateral','ipsilateral')
annotation('textbox', [.055, (1/3)-0.1+0.01, .15, .05], 'string', '5 - 16 kHz','fontsize',18,'LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle')

% save figure
figlen = 8;
width = 4*figlen;
height = 7*figlen;
set(gcf,'Units','centimeters','PaperPosition',[0 0 width height],'PaperSize',[width height]);
saveas(gcf,'data/hmdpert_output/sd_difference_interpolated.png')  


%% CREATE INVERSE FILTERS
figure('Name','correction filters','NumberTitle','off','WindowStyle','docked');
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

    plot(f,20*log10(iH))
    
    set(gca,'xscale','log')
    grid on
    xlim([10 24000])
    ylim([-12 12])
    
    subplot(2,1,2)
    hold on
    plot(invh)
    
    model_interp(i).invh = invh;
end

save('data/hmdpert_output/model_interp.mat', 'model_interp')


%% functions
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

function plotAzElM(az,el,val,lim,clbl,texta,textb)
    step = 0.5; % 0.5
    azimuth = -180:step:180;
    elevation = -90:step:90;
    azimuth = azimuth + randn(size(azimuth))*0.01;
    elevation = elevation + randn(size(elevation))*0.01;
    azimuth = min(max(azimuth,-180),180);
    elevation = min(max(elevation,-90),90);

    for i = 1:length(azimuth)
        for j = 1:length(elevation)
            dist = distance(elevation(j),azimuth(i),el,az);
            idx = find(dist <= min(dist) + 0.001);
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

    % labels
    fsize = 16;
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
    
    text(-2,1.25,texta,'color',fcolor,'fontsize',fsize,'rotation',0,'horizontalalignment','center','verticalalignment','middle');
    text(2,1.25,textb,'color',fcolor,'fontsize',fsize,'rotation',0,'horizontalalignment','center','verticalalignment','middle');
end
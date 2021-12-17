close all
clear

addpath('tools/')
addpath('tools/invFIR/')
addpath('tools/VoronoiSphere/')
addpath('../API_MO/API_MO/')

subjectdir = 'data/20210923-RIG-KEMAR-NOHMD/';

% subjectdir = 'data/20201217-122pt-2.5m-dayton_vt/';
% subjectdir = 'data/20201217-122pt-2.5m-canford_vt/';
% subjectdir = 'data/20211012-q2_tr/';
% subjectdir = 'data/20211105-A-Jan/';
% subjectdir = 'data/20211126-XR-TR/';
% subjectdir = 'data/20211126-Gavin/';

load([subjectdir 'irBank.mat'])
mkdir([subjectdir 'figures/'])
plotMagnitudes(irBank, '1-measured', [subjectdir 'figures/'])

% time domain windowing
plotting = 'false';
irBank = winIRs(irBank, plotting, [subjectdir 'figures/windowing/']); % set 'true' to save plots
plotMagnitudes(irBank, '2-win', [subjectdir 'figures/'])

% calculate FF measurement inv filter and normalize all hrirs
% do the low-frequency extension
irBank = normalizeIRs(irBank, 'true');
plotMagnitudes(irBank, '3-raw', [subjectdir 'figures/'])

% diffuse-field equalization
dfe_enabled = true;
irBank = dfeHRIRs(irBank, dfe_enabled, 'true', [subjectdir 'figures/']);
plotMagnitudes(irBank, '4-dfe', [subjectdir 'figures/'])

% save ambix config file
saveAsAmbix(irBank, subjectdir)

% save sofa file
saveAsSofa(irBank, subjectdir)

function IRbank = winIRs(IRbank, plotting, save_fig_folder)
    % define window
    win = hann(300).^4;
    win = [win(1:end/2); ones(80,1); win(end/2+1:end);];
    winshift = 160; % how many samples the window should be shifted forward from the peak
    Nwin = length(win);
    if mod(Nwin,2) ~= 0
        disp('wrong window length, must be even')
    end
    
    for i = 1:length(IRbank)
        irLeft = [zeros(Nwin,1); IRbank(i).fullIR(:,1)];
        irRight = [zeros(Nwin,1); IRbank(i).fullIR(:,2)];
        Fs = IRbank(i).Fs;
        [IRbank(i).ITD, maxL, maxR] = getITD(irLeft,irRight,Fs);

        % apply window
        winstart = maxL - winshift;
        winend = winstart + Nwin - 1;
        irLeft(1:winstart-1) = 0;
        irLeft(winstart:winend) = irLeft(winstart:winend) .* win;
        irLeft(winend+1:end) = 0;
        
        winstart = maxR - winshift;
        winend = winstart + Nwin - 1;
        irRight(1:winstart-1) = 0;
        irRight(winstart:winend) = irRight(winstart:winend) .* win;
        irRight(winend+1:end) = 0;
        
        % mean time of arrival for both ear signals expressed in samples
        toasmp = round(mean([maxL, maxR]));
        
        % cut irs
        preSamples = floor(0.002 * Fs);
        afterSamples = Nwin-preSamples;
        IRbank(i).winIR(:,1) = irLeft(toasmp-preSamples+1:toasmp+afterSamples);
        IRbank(i).winIR(:,2) = irRight(toasmp-preSamples+1:toasmp+afterSamples);
        IRbank(i).maxL = maxL - Nwin; % because Nwin of zeroes was added at the beginning
        IRbank(i).maxR = maxR - Nwin;
        IRbank(i).toasmp = toasmp - Nwin;
        IRbank(i).winMaxL = maxL - toasmp + preSamples;
        IRbank(i).winMaxR = maxR - toasmp + preSamples;
    end
    
    % plot ITD
    for i = 1:length(IRbank)
        IRbank(i).ITDwin = (IRbank(i).maxR - IRbank(i).maxL)  * 10^6 / IRbank(i).Fs;
    end
    figure('Name','ITD','NumberTitle','off','WindowStyle','docked');
    tiledlayout(1,2)
    lim = [-1000 1000];
    nexttile
    hold on
    plotAzEl([IRbank.azimuth],[IRbank.elevation],[IRbank.ITD], lim)
    nexttile
    hold on
    plotAzEl([IRbank.azimuth],[IRbank.elevation],[IRbank.ITDwin], lim)
    
    % plot distance
    figure('Name','distance','NumberTitle','off','WindowStyle','docked');
    scatter([IRbank.distance],[IRbank.toasmp] * 343/Fs)
    xlabel('measured distance by HMD (m)')
    ylabel('measured distance acoustically (m)')
    
    % correct gains
    ref_dist = 1.5;
    for i = 1:length(IRbank)
        if ~isnan(IRbank(i).distance)
            gain_lin = IRbank(i).distance/ref_dist;
            IRbank(i).gain = 20*log(gain_lin);
            IRbank(i).winIR = IRbank(i).winIR * gain_lin;
        end
    end
    
    if strcmp(plotting,'true')
        figure('Name','ETC + win','NumberTitle','off','WindowStyle','docked');
        range = [-20 20];
        mkdir(save_fig_folder)
        for i = 1:length(IRbank) % this will produce lots of graphs
%         for i = randperm(length(IRbank),20) % 20 random directions
            subplot(2,2,1)
            hold on
            yyaxis left
            cla
            plot(IRbank(i).fullIR(:,1),'-b')
            xline(IRbank(i).maxL, '--k')
            xline(IRbank(i).toasmp, '--g')
            xlabel('Samples')
            ylabel('Amplitude')
            ylim(range)
            yyaxis right
            cla
            plot(IRbank(i).maxL-winshift:IRbank(i).maxL-winshift+Nwin-1,win,'--')
            ylabel('Linear gain')
            xlim([IRbank(i).maxL-Nwin IRbank(i).maxL+Nwin])
            title(['Left Raw HRIR' ' azi ' num2str(IRbank(i).azimuth) ' ele ' num2str(IRbank(i).elevation)])

            subplot(2,2,2)
            cla
            hold on
            plot(IRbank(i).winIR(:,1),'-b')
            xlabel('Samples')
            ylabel('Amplitude')
            ylim(range);
            xlim([0 size(IRbank(i).winIR,1)])
            title('Left Windowed HRIR')
            
            subplot(2,2,3)
            hold on
            yyaxis left
            cla
            plot(IRbank(i).fullIR(:,2),'-b')
            xline(IRbank(i).maxR, '--k')
            xline(IRbank(i).toasmp, '--g')
            xlabel('Samples')
            ylabel('Amplitude')
            ylim(range)
            yyaxis right
            cla
            plot(IRbank(i).maxR-winshift:IRbank(i).maxR-winshift+Nwin-1,win,'--')
            ylabel('Linear gain')
            xlim([IRbank(i).maxR-Nwin IRbank(i).maxR+Nwin])
            title('Right Raw HRIR')

            subplot(2,2,4)
            cla
            hold on
            plot(IRbank(i).winIR(:,2),'-b')
            xlabel('Samples')
            ylabel('Amplitude')
            ylim(range);
            xlim([0 size(IRbank(i).winIR,1)])
            title('Right Windowed HRIR')
            
            % save figure
            figlen = 8;
            width = 4*figlen;
            height = 3*figlen;
            set(gcf,'Units','centimeters','PaperPosition',[0 0 width height],'PaperSize',[width height]);
            saveas(gcf,[save_fig_folder 'windowing_fig' num2str(i)  '.jpg'])
%             close
        end
    end
end

function irBank = normalizeIRs(irBank, plotting)
    % find the free-field measurements and calculate inverse filters
    for i = find([irBank.ref])
        Fs = irBank(i).Fs;
        FFIR(:,1) = flattenMagHF(irBank(i).winIR(:,1), Fs);
        FFIR(:,2) = flattenMagHF(irBank(i).winIR(:,2), Fs);
        Nfft = length(irBank(i).winIR);
        irBank(i).invh = invFIR('minphase', FFIR, Nfft, 0, Nfft, [0 Fs/2], [0 0], 1, Fs);
    end
    
    if strcmp(plotting,'true')
        % plot freefield and inverse filter
        plots_num = length(find([irBank.ref]));
        figure('Name','P0 + inv resp','NumberTitle','off','WindowStyle','docked');
        tiledlayout(ceil(sqrt(plots_num)),floor(sqrt(plots_num)))
        for i = find([irBank.ref])
            nexttile
            hold on
            box on

            [f,mag] = getMagnitude(irBank(i).winIR(:,1),irBank(i).Fs,'log');
            plot(f,mag,'-g','LineWidth',1);

            [f,mag] = getMagnitude(irBank(i).winIR(:,2),irBank(i).Fs,'log');
            plot(f,mag,'-r','LineWidth',1);

            [f,mag] = getMagnitude(irBank(i).invh(:,1),irBank(i).Fs,'log');
            plot(f,mag,'--g','LineWidth',2);

            [f,mag] = getMagnitude(irBank(i).invh(:,2),irBank(i).Fs,'log');
            plot(f,mag,'--r','LineWidth',2);

            set(gca,'xscale','log')
            grid on
            xlim([10 Fs/2]);
            % ylim([-80 -30]);
%             legend('Left channel', 'Right channel','','', 'Left inverse filter', 'Right inverse filter','location','northwest')
            xlabel('Frequency (Hz)');
            ylabel('Magnitude (dB)');
        end

%         % plot freefield impulse response
%         figure('Name','FF IR','NumberTitle','off','WindowStyle','docked');
%         subplot(2,1,1)
%         hold on
%         box on
%         plot(FFIR(:,1),'g')
%         plot(FFIR(:,2),'r')
%         
%         % plot inverese filter for freefield impulse response
% %         figure('Name','FF inv filter IR','NumberTitle','off','WindowStyle','docked');
%         subplot(2,1,2)
%         hold on
%         box on
%         plot(invh(:,1),'g')
%         plot(invh(:,2),'r')
    end
    
    %% filter all windowed hrirs
    for i = 1:length(irBank)
        % find the reference measurement
        ref_idx = [irBank.ref] == 1 & [irBank.lspk] == irBank(i).lspk;
        irBank(i).rawHRIR(:,1) = conv(irBank(ref_idx).invh(:,1),irBank(i).winIR(:,1));
        irBank(i).rawHRIR(:,2) = conv(irBank(ref_idx).invh(:,2),irBank(i).winIR(:,2));  
    end
    
    %% add LF extension
    for i = 1:length(irBank)
        irBank(i).rawHRIR(:,1) = LFextension(irBank(i).rawHRIR(:,1), irBank(i).winMaxL, Fs);
        irBank(i).rawHRIR(:,2) = LFextension(irBank(i).rawHRIR(:,2), irBank(i).winMaxR, Fs);  
    end
    
    %% cut equalized IRs
    figure('Name','raw IR truncation','NumberTitle','off','WindowStyle','docked');
    hold on
    cut_start = 1;
    cut_end = cut_start + 255 + 256;
    for i = 1:length(irBank)
        plot(irBank(i).rawHRIR(:,1),'b')
        plot(irBank(i).rawHRIR(:,2),'r')
    end
    xline(cut_start,'--k')
    xline(cut_end,'--k')
    
    for i = 1:length(irBank)
        irBank(i).rawHRIR = irBank(i).rawHRIR(cut_start:cut_end,:);
    end
    
end

function irBank = dfeHRIRs(irBank, dfe_enabled, plotting, save_fig_folder)
    Fs = unique([irBank.Fs]);
    
    % remove the ff measurement
%     IRbank(isnan([IRbank.azimuth])) = [];
    irBank([irBank.ref]) = [];

    % calculate weights based on solid angle of each cell
    azel = [[irBank.azimuth]',[irBank.elevation]'];
    azel = azel + randn(size(azel))*0.001; % add some noise to smooth out plotting
    [xyz(1,:), xyz(2,:), xyz(3,:)] = sph2cart(deg2rad(azel(:,1)),deg2rad(azel(:,2)),1);
    xyz = xyz ./ sqrt(sum(xyz.^2,1));
    [~, ~, voronoiboundary, s] = voronoisphere(xyz);
    
    sn = s./sum(s);
    
    % calculate average magnitude for left & right
%     Nfft = 4096;    
    for i = 1:length(irBank)
        if i == 1
            Nfft = size(irBank(i).rawHRIR,1);
            mag_ir_avgL = zeros(Nfft,1);
            mag_ir_avgR = zeros(Nfft,1);
        end
        mag_ir_avgL = mag_ir_avgL + (abs(fft(irBank(i).rawHRIR(:,1), Nfft)).^2) * sn(i);
        mag_ir_avgR = mag_ir_avgR + (abs(fft(irBank(i).rawHRIR(:,2), Nfft)).^2) * sn(i);
    end
    mag_ir_avgL = sqrt(mag_ir_avgL);
    mag_ir_avgR = sqrt(mag_ir_avgR);

    % back to time domain
    ir_avgL = ifft(mag_ir_avgL,'symmetric');
    ir_avgR = ifft(mag_ir_avgR,'symmetric');
    ir_avgL=circshift(ir_avgL,Nfft/2,1);
    ir_avgR=circshift(ir_avgR,Nfft/2,1);
    
    % flatten
    ir_avgL = flattenMagHF(ir_avgL, Fs);
    ir_avgR = flattenMagHF(ir_avgR, Fs);
    
    % inv fir Nfft and filter length
    Nfft = length(ir_avgL);
    if dfe_enabled   
        % create dfe filters
        dfeL = invFIR('minphase', ir_avgL, Nfft, 0, Nfft, [0 Fs/2], [0 0], 1, Fs);
        dfeR = invFIR('minphase', ir_avgR, Nfft, 0, Nfft, [0 Fs/2], [0 0], 1, Fs);
    else
        % use "flat" filters
        dfeL = [1; zeros(Nfft-1,1)];
        dfeR = [1; zeros(Nfft-1,1)];
    end

    %% filter hrirs with dfe filters
    for i = 1:length(irBank)
        irBank(i).dfeHRIR(:,1) = conv(dfeL,irBank(i).rawHRIR(:,1));
        irBank(i).dfeHRIR(:,2) = conv(dfeR,irBank(i).rawHRIR(:,2));  
    end
    
    %% cut equalized IRs
    figure('Name','dfe IR truncation','NumberTitle','off','WindowStyle','docked');
    hold on
    cut_start = 1;
    cut_end = cut_start + 255;
    for i = 1:length(irBank)
        plot(irBank(i).dfeHRIR(:,1),'b')
        plot(irBank(i).dfeHRIR(:,2),'r')
    end
    xline(cut_start,'--k')
    xline(cut_end,'--k')
    
    for i = 1:length(irBank)
        irBank(i).dfeHRIR = irBank(i).dfeHRIR(cut_start:cut_end,:);
    end
    
    if plotting == "true"        
        
        figure('Name','voronoi diagram','NumberTitle','off','WindowStyle','docked');
        hold on
        n = size(xyz,2);
        w = s - min(s);
        w = w ./ max(w);
        w = round(w*255)+1;
            plot3(xyz(1,:),xyz(2,:),xyz(3,:),'ko');
%         text(xyz(1,:),xyz(2,:),xyz(3,:),num2str(w))
        clmap = cool();
        ncl = size(clmap,1);
        for k = 1:n
            X = voronoiboundary{k};
            cl = clmap(w(k),:);
            fill3(X(1,:),X(2,:),X(3,:),cl,'EdgeColor','w');
        end
        axis('equal');
        view(40,10)
        axis([-1 1 -1 1 -1 1]);
        
        % plot avereged hrirs and dfe filter response 
        figure('Name','avg resp + dfe filter','NumberTitle','off','WindowStyle','docked');
        hold on
        box on
        grid on
        Nfft = 1024;
        f = (Fs/Nfft:Fs/Nfft:Fs)';
        mag = 20*log10(abs(fft(ir_avgL, Nfft)));
        plot(f,mag,'-g','LineWidth',1);
        mag = 20*log10(abs(fft(ir_avgR, Nfft)));
        plot(f,mag,'-r','LineWidth',1);
        mag = 20*log10(abs(fft(dfeL, Nfft)));
        plot(f,mag,'--g','LineWidth',2);
        mag = 20*log10(abs(fft(dfeR, Nfft)));
        plot(f,mag,'--r','LineWidth',2);
                
        set(gca,'xscale','log')
        xlim([10 Fs/2]);
        ylim([-20 20]);
        ylabel('Magnitude (dB)')
        xlabel('Frequency (Hz)')
        legend('Left ear avg', 'Rigth ear avg', 'Left DFE filter', 'Right DFE filter')
        title('DFE filter')
        
        % save figure
        figlen = 8;
        width = 4*figlen;
        height = 3*figlen;
        set(gcf,'Units','centimeters','PaperPosition',[0 0 width height],'PaperSize',[width height]);
        mkdir(save_fig_folder)
        saveas(gcf,[save_fig_folder 'dfe.jpg'])        
%         close
        
%         % equalized magnitudes
%         figure('Name','Magnitudes','NumberTitle','off','WindowStyle','docked');
%         hold on
%         for i = find([IRbank.azimuth] == 0 & [IRbank.elevation] == 0) %1:length(IRbank)
%             Nfft = 1024;
%             Fs = IRbank(i).Fs;
%             f = (Fs/Nfft:Fs/Nfft:Fs)';
% %             raw hrir
%             plot(f,20*log10(abs(fft(IRbank(i).rawhrir(:,1),Nfft))),'-','linewidth',0.5)
%             plot(f,20*log10(abs(fft(IRbank(i).rawhrir(:,2),Nfft))),'-','linewidth',0.5)
% 
% %             windowed hrir
%             plot(f,20*log10(abs(fft(IRbank(i).winhrir(:,1),Nfft))),'-','linewidth',0.5)
%             plot(f,20*log10(abs(fft(IRbank(i).winhrir(:,2),Nfft))),'-','linewidth',0.5)
% 
%             % dfe hrir
%             plot(f,20*log10(abs(fft(IRbank(i).dfehrir(:,1),Nfft))),'-','linewidth',0.5)
%             plot(f,20*log10(abs(fft(IRbank(i).dfehrir(:,2),Nfft))),'-','linewidth',0.5)
%         end
%         box on
%         set(gca,'xscale','log')
%         xticks([10 100 1000 10000])
%         xlim([40 Fs/2])
%         ylim([-60 20])
%         xlabel('Frequency (Hz)')
%         ylabel('Magnitude (dB)')
    end
end

function ir_fixed = flattenMagHF(ir, Fs)
    Nfft = length(ir);
    f = (Fs/Nfft:Fs/Nfft:Fs);
    mag_target = abs(fft(ir, Nfft));
    
    % establish amplitude level
    idmin = find(f >= 16000, 1 );           % f min
    idmax = find(f <= 20000, 1, 'last');    % f max
    hfext_amp = mean(mag_target(idmin:idmax));
    
    % define "crossover" region
    idmin = find(f >= 18000, 1 );           % f min
    idmax = find(f <= 20000, 1, 'last');    % f max
    
    % flatten above fmax
    win = hann(2*(idmax-idmin+1));
    win1 = win(length(win)/2+1:end);
    win2 = win(1:length(win)/2);
    mag_fixed = mag_target;
    mag_fixed(idmin:idmax) = mag_target(idmin:idmax) .* win1 + hfext_amp .* win2;
    mag_fixed(idmax+1:end) = hfext_amp;
    
    % back to time domain
    ir_fixed = ifft(mag_fixed,'symmetric');
    ir_fixed=circshift(ir_fixed,Nfft/2,1);
%     ir_fixed=0.5*(1-cos(2*pi*(1:Nfft)'/(Nfft+1))).*ir_fixed;
end

function lfe_ir = LFextension(ir, shift,  Fs)
    lfe_ir = zeros(length(ir),1);
    lfe_ir(shift) = 1;
    xfreq = 250;
    filter_order = 2;    
    [B_highpass, A_highpass] = butter( filter_order, xfreq/Fs*2, 'high' );            
    [B_lowpass,  A_lowpass ] = butter( filter_order, xfreq/Fs*2, 'low'  );
    output_low = filter(B_lowpass, A_lowpass, lfe_ir);
    output_low = filter(B_lowpass, A_lowpass, output_low);
    output_high = filter(B_highpass, A_highpass, ir);
    output_high = filter(B_highpass, A_highpass, output_high);
    
    lfe_ir = output_low + output_high;
end

function plotMagnitudes(irBank, type, save_fig_folder)
    plots_num = 50;
    fig_num = ceil(length(irBank)/plots_num);

    for i = 1:fig_num
        %% plot full ir spectrum
        figure('Name',[type 'mag-' num2str(i)],'NumberTitle','off','WindowStyle','docked');
        tiledlayout(ceil(sqrt(plots_num)),floor(sqrt(plots_num)))
        for j = (i-1)*plots_num+1:min(i*plots_num,length(irBank))
            Fs = irBank(j).Fs;
            nexttile
            hold on
            if strcmp(type,'1-measured')
                [f,mag] = getMagnitude(irBank(j).fullIR(:,1),Fs,'log');
                plot(f,mag,'-b','LineWidth',1);
                [f,mag] = getMagnitude(irBank(j).fullIR(:,2),Fs,'log');
                plot(f,mag,'-r','LineWidth',1);
            elseif strcmp(type,'2-win')
                [f,mag] = getMagnitude(irBank(j).winIR(:,1),Fs,'log');
                plot(f,mag,'-b','LineWidth',1);
                [f,mag] = getMagnitude(irBank(j).winIR(:,2),Fs,'log');
                plot(f,mag,'-r','LineWidth',1);
            elseif strcmp(type,'3-raw')
                [f,mag] = getMagnitude(irBank(j).rawHRIR(:,1),Fs,'log');
                plot(f,mag,'-b','LineWidth',1);
                [f,mag] = getMagnitude(irBank(j).rawHRIR(:,2),Fs,'log');
                plot(f,mag,'-r','LineWidth',1);
            elseif strcmp(type,'4-dfe')
                [f,mag] = getMagnitude(irBank(j).dfeHRIR(:,1),Fs,'log');
                plot(f,mag,'-b','LineWidth',1);
                [f,mag] = getMagnitude(irBank(j).dfeHRIR(:,2),Fs,'log');
                plot(f,mag,'-r','LineWidth',1);
            end

    %         legend('Left Ear','Right Ear','location','southwest')
    %         xlim([100 20000]);
            xlim([20 Fs/2]);
            ylim([-40 40]);
            xlabel('Frequency (Hz)');
            ylabel('Magnitude (dB)');
            title(['az: ' num2str(irBank(j).azimuth) ', el: ' num2str(irBank(j).elevation) ', dist: ' num2str(irBank(j).distance)]);
            set(gca,'xscale','log')
        end

        % save figure
        figlen = 16;
        width = 4*figlen;
        height = 3*figlen;
        set(gcf,'PaperPosition',[0 0 width height],'PaperSize',[width height]);
        saveas(gcf,[save_fig_folder 'measured_mag_' type '-' num2str(i) '.jpg'])
%         close
    end

end

function saveAsAmbix(IRbank, subjectdir)
    addpath('tools/')
    addpath('../adt/')
    adt_initialize

    % remove the ff measurement
    IRbank(isnan([IRbank.azimuth])) = [];
    
    % save wav files
    mkdir(subjectdir,'ambix_wav') 
    for i = 1:length(IRbank)
        filename = ['azi_' num2str(IRbank(i).azimuth,'%.2f') '_ele_' num2str(IRbank(i).elevation,'%.2f') '.wav'];
        disp(filename);
        IRbank(i).filename = filename;
        gain = 0.5; %10^(-3/20);
%         audiowrite([subjectdir '/ambix_wav/' filename], gain*[IRbank(i).hrirLeft IRbank(i).hrirRight], IRbank(i).Fs)
        audiowrite([subjectdir '/ambix_wav/' filename], gain * IRbank(i).dfeHRIR, IRbank(i).Fs)
    end
    
    % create decoder presets
    % preset 1
    layout_list(1).order = 3;
    layout_list(1).layout = '26Leb';
    ls = getLebedevSphere(26);
    dirs = [];
    [dirs(:,1), dirs(:,2), ~] = cart2sph(ls.x,ls.y,ls.z);
    layout_list(1).dirs = rad2deg(dirs);
    
    % preset 2
    layout_list(2).order = 5;
    layout_list(2).layout = '50Leb';
    ls = getLebedevSphere(50);
    dirs = [];
    [dirs(:,1), dirs(:,2), ~] = cart2sph(ls.x,ls.y,ls.z);
    layout_list(2).dirs = rad2deg(dirs);
    
    for i = 1:length(layout_list)
        order = layout_list(i).order;
        layout = layout_list(i).layout;
        dirs = layout_list(i).dirs;
        ambi = struct;
        for j = 1:length(dirs)
           % fine nearest measured point
           [~,idx] = min(distance([IRbank.elevation],-[IRbank.azimuth],dirs(j,2),dirs(j,1)));
           ambi(j).az = -IRbank(idx).azimuth;
           ambi(j).el = IRbank(idx).elevation;
           ambi(j).dist = IRbank(idx).distance;
           ambi(j).filename = IRbank(idx).filename;
        end

        S = SPKR_ARRAY([ambi.az],[ambi.el],[ambi.dist]);
        mkdir([subjectdir 'ambi_dec/'])
        out_path = [subjectdir 'ambi_dec/xr-hrtf-' num2str(order) 'OA-' layout];
        do_plots = false;
        decoder_type = 2;
        [D, S, M, C] = ambi_run_pinv(S,order,[],out_path,do_plots,[],0,decoder_type);
        mtx = M.lf; % get the basic matrix
        
        % write the ambix config file
        out_path = [subjectdir 'ambix_wav/xr-hrtf-' num2str(order) 'OA-' layout];  
        fid = fopen([out_path '.config'],'w');
        fprintf(fid, '// Ambix config file\n');
        fprintf(fid, [...
            '\n', ...
            '#GLOBAL\n', ...
            '/debug_msg %s\n', ...
            '/coeff_scale %s\n', ...
            '/coeff_seq  %s\n', ...
            '/flip %i\n', ...  % 1 negates y axis
            '/dec_mat_gain %f\n', ... % set to 1 per MK
            '#END\n'], ...
            D.description, lower(D.coeff_scale), lower(D.coeff_order), ...
            0, ...  % flip
            1  ...  % dec_mat_gain
            );
        
        fprintf(fid,'\n#HRTF\n');
        for j = 1:length(ambi)
            fprintf(fid,[ambi(j).filename '\n']);
        end
        fprintf(fid,'#END\n');
        fprintf(fid, '\n#DECODERMATRIX\n');
        for j = 1:size(mtx,1)
            fprintf(fid, '\t% f', mtx(j,:));
            fprintf(fid, '\n');
        end
        fprintf(fid, '#END\n');
        fclose(fid);
    end
end

function [ val ] = SPKR_ARRAY(azis,eles,dist)
    val.name = 'xr-hrtf';
    % X Y Z
    [S(:,1),S(:,2),S(:,3)] = sph2cart(deg2rad(azis),deg2rad(eles),dist);    
    % spherical coordinate representation
    [val.az, val.el, val.r] = cart2sph(S(:,1),S(:,2),S(:,3));
    
    % unit vectors
    [val.x, val.y, val.z] = sph2cart(val.az, val.el, 1);
   
    for i = 1:length(azis)
        val.id(i) = {['spk' num2str(i)]};
    end
end

function saveAsSofa(IRbank, subjectdir)
    % Start SOFA
    SOFAstart
    
    % Get an empy conventions structure
    Obj = SOFAgetConventions('SimpleFreeFieldHRIR');
    
    % remove the ff measurement
    IRbank(isnan([IRbank.azimuth])) = [];

    for i = 1:length(IRbank)
%        hrirs(i,1,:) = IRbank(i).hrirLeft;
%        hrirs(i,2,:) = IRbank(i).hrirRight;
       hrirs(i,:,:) = IRbank(i).dfeHRIR';
       azi = IRbank(i).azimuth;
       ele = IRbank(i).elevation;
       dist = IRbank(i).distance;

       if(azi < 0)
           azi = azi + 360;
       end
       azi = azi * -1;

       Obj.SourcePosition(i,:)=[azi ele dist];
    end
    
    Obj.Data.IR = hrirs;
    Obj.Data.SamplingRate = IRbank(1).Fs;
    
    % Update dimensions
    Obj=SOFAupdateDimensions(Obj);

    % %% Fill with attributes
    % Obj.GLOBAL_ListenerShortName = 'KEMAR';
    % Obj.GLOBAL_History = 'created with a script';
    % Obj.GLOBAL_DatabaseName = 'none';
    % Obj.GLOBAL_ApplicationName = 'Demo of the SOFA API';
    % Obj.GLOBAL_ApplicationVersion = SOFAgetVersion('API');
    % Obj.GLOBAL_Organization = 'Acoustics Research Institute';
    % Obj.GLOBAL_AuthorContact = 'piotr@majdak.com';
    % Obj.GLOBAL_Comment = 'Contains simple pulses for all directions';

    %% save the SOFA file
    % Data compression (0..uncompressed, 9..most compressed)
    compression=1; % results in a nice compression within a reasonable processing time
    SOFAfn=fullfile(subjectdir,'xr-hrtf.sofa');
    disp(['Saving:  ' SOFAfn]);
    SOFAsave(SOFAfn, Obj, compression);
end

function plotAzEl(az,el,val, lim)
    azimuth = -180:1:180;
    elevation = -90:1:90;
    
    % get absolute val
%     val = abs(val);

    for i = 1:length(azimuth)
        for j = 1:length(elevation)
            dist = distance(elevation(j),azimuth(i),el,az);
            idx = find(dist == min(dist));
            if length(idx) > 1
                value(j,i) = median(val(idx));
            else
                value(j,i) = val(idx);
            end

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
%     title('MFSD')
    colorbar
end
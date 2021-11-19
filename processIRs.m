close all
clear

addpath('tools/invFIR/')
addpath('tools/VoronoiSphere/')
addpath('../API_MO/API_MO/')

% subjectdir = 'data/20201217-122pt-2.5m-dayton_vt/';
% subjectdir = 'data/20201217-122pt-2.5m-canford_vt/';
% subjectdir = 'data/20211012-q2_tr/';
subjectdir = 'data/20211104-q2_tr/';
% subjectdir = 'data/20211105-A/';

load([subjectdir 'irBank.mat'])
plotMagnitudes(irBank, '1-measured', [subjectdir 'figures/'])

% time domain windowing
irBank = winIRs(irBank, 'true', [subjectdir 'figures/windowing/']); % set 'true' to save plots
plotMagnitudes(irBank, '2-win', [subjectdir 'figures/'])

% calculate FF measurement inv filter and normalize all hrirs
irBank = normalizeIRs(irBank, 'true');
plotMagnitudes(irBank, '3-raw', [subjectdir 'figures/'])

% do the low-frequency extension

% diffuse-field equalization
dfe_enabled = true;
irBank = dfeHRIRs(irBank, dfe_enabled, 'true', [subjectdir 'figures/']);
plotMagnitudes(irBank, '4-dfe', [subjectdir 'figures/'])

% save ambix config file
saveAsAmbix(irBank, subjectdir)

% save sofa file
saveAsSofa(irBank, subjectdir)

function IRbank = winIRs(IRbank, plotting, save_fig_folder)
%     % default option
%     Nwin = 512;
%     win = hann(Nwin).^6;

%     % new option
%     win = hann(200).^4;
%     win = [win(1:end/2); ones(56,1); win(end/2+1:end);];
%     Nwin = length(win);

    % new option 2
    win = hann(300).^4;
    win = [win(1:end/2); ones(100,1); win(end/2+1:end);];
    Nwin = length(win);
    
    for i = 1:length(IRbank)
        irLeft = [zeros(Nwin,1); IRbank(i).fullIrLeft; zeros(Nwin,1)];
        irRight = [zeros(Nwin,1); IRbank(i).fullIrRight; zeros(Nwin,1)];
        Fs = IRbank(i).Fs;
        [~,maxL] = max(abs(irLeft));
        [~,maxR] = max(abs(irRight));
        
        % check if maxL and maxR are spaced from each other by more than
        % 900 us
        n = floor(0.0009 * Fs); %43 samples@48k, approximate maximum itd
        while abs(maxL - maxR) > n
            if maxL - maxR > n
                [~,maxL] = max(abs(irLeft(1:maxR+n)));
            elseif maxR - maxL > n
                [~,maxR] = max(abs(irRight(1:maxL+n)));
            end
        end
        
        % create empty vectors
        IRbank(i).winIrLeft = zeros(2*Nwin,1);
        IRbank(i).winIrRight = zeros(2*Nwin,1);
        
        % calculate the lag of hrirs in reference to the mean TOA
        shiftL = maxL - floor(mean([maxL maxR]));
        shiftR = maxR - floor(mean([maxL maxR]));
        
        % do the windowing and copy output to the new vectors with
        % appropriate sample index shifts
        IRbank(i).winIrLeft(Nwin-Nwin/2+1+shiftL:Nwin+Nwin/2+shiftL) = irLeft(maxL-Nwin/2+1:maxL+Nwin/2) .* win;
        IRbank(i).winIrRight(Nwin-Nwin/2+1+shiftR:Nwin+Nwin/2+shiftR) = irRight(maxR-Nwin/2+1:maxR+Nwin/2) .* win;
        preSamples = floor(0.001 * Fs); %48 samples @ 48kHz
        IRbank(i).winIrLeft = IRbank(i).winIrLeft(end/2-preSamples+1:end-preSamples);
        IRbank(i).winIrRight = IRbank(i).winIrRight(end/2-preSamples+1:end-preSamples);
        IRbank(i).maxL = maxL - Nwin; % because Nwin of zeroes was added at the beginning and end
        IRbank(i).maxR = maxR - Nwin;
    end
    
    if strcmp(plotting,'true')
        figure('Name','ETC + win','NumberTitle','off','WindowStyle','docked');
%         figure
        range = [-20 20];
        mkdir(save_fig_folder)
        for i = 1:length(IRbank) % this will produce lots of graphs
%         for i = randperm(length(IRbank),20) % 20 random directions
            subplot(2,2,1)
            hold on
            yyaxis left
            cla
%             plot(20*log10(sqrt(IRbank(i).fullIrLeft.^2)),'-b')
            plot(IRbank(i).fullIrLeft,'-b')
            xline(IRbank(i).maxL, '--k')
            xlabel('Samples')
            ylabel('Amplitude')
            ylim(range)
            yyaxis right
            cla
            plot(IRbank(i).maxL-Nwin/2+1:IRbank(i).maxL+Nwin/2,win,'--')
            ylabel('Linear gain')
            xlim([IRbank(i).maxL-Nwin IRbank(i).maxL+Nwin])
            title(['Left Raw HRIR' ' azi ' num2str(IRbank(i).azimuth) ' ele ' num2str(IRbank(i).elevation)])

            subplot(2,2,2)
            cla
            hold on
            plot(IRbank(i).winIrLeft,'-b')
            xlabel('Samples')
            ylabel('Amplitude')
            ylim(range);
            xlim([0 Nwin])
            title('Left Windowed HRIR')
            
            subplot(2,2,3)
            hold on
            yyaxis left
            cla
            plot(IRbank(i).fullIrRight,'-b')
            xline(IRbank(i).maxR, '--k')
            xlabel('Samples')
            ylabel('Amplitude')
            ylim(range)
            yyaxis right
            cla
            plot(IRbank(i).maxR-Nwin/2+1:IRbank(i).maxR+Nwin/2,win,'--')
            ylabel('Linear gain')
            xlim([IRbank(i).maxR-Nwin IRbank(i).maxR+Nwin])
            title('Right Raw HRIR')

            subplot(2,2,4)
            cla
            hold on
            plot(IRbank(i).winIrRight,'-b')
            xlabel('Samples')
            ylabel('Amplitude')
            ylim(range);
            xlim([0 Nwin])
            title('Right Windowed HRIR')
            
            % save figure
            figlen = 8;
            width = 4*figlen;
            height = 3*figlen;
            set(gcf,'Units','centimeters','PaperPosition',[0 0 width height],'PaperSize',[width height]);
            saveas(gcf,[save_fig_folder 'windowing_fig' num2str(i)  '.jpg'])        
        end
        close
    end
end

function irBank = normalizeIRs(irBank, plotting)
    % pick the free-field measurement and calculate inverse filter
    i = length(irBank);
    Fs = irBank(i).Fs;
    Nfft = length(irBank(i).winIrLeft)*8;
    FFIR = [irBank(i).winIrLeft irBank(i).winIrRight];
    invh = invFIR('minphase', FFIR, Nfft, 0, Nfft, [10 Fs/2], [120 0], 1, Fs);
    
    % low pass the inverse filter
    lpFilt = designfilt('lowpassiir','FilterOrder',2, ...
         'PassbandFrequency',22000,'PassbandRipple',0.5, ...
         'SampleRate',Fs);
    invh = filter(lpFilt,invh);
    
    if strcmp(plotting,'true')
        % plot freefield and inverse filter
        figure('Name','P0 + inv resp','NumberTitle','off','WindowStyle','docked');
        hold on
        box on

        fullIR = irBank(i).winIrLeft;
        Nfft = length(fullIR);
        f = (Fs/Nfft:Fs/Nfft:Fs)';
        mag = 20*log10(abs(fft(fullIR, Nfft)));
        plot(f,mag,'-g','LineWidth',1);
        
        fullIR = invh(:,1);
        Nfft = length(fullIR);
        f = (Fs/Nfft:Fs/Nfft:Fs)';
        mag = 20*log10(abs(fft(fullIR, Nfft)));
        plot(f,mag,'--g','LineWidth',2);

        fullIR = irBank(i).winIrRight;
        Nfft = length(fullIR);
        f = (Fs/Nfft:Fs/Nfft:Fs)';
        mag = 20*log10(abs(fft(fullIR, Nfft)));
        plot(f,mag,'-r','LineWidth',1);
        
        fullIR = invh(:,2);
        Nfft = length(fullIR);
        f = (Fs/Nfft:Fs/Nfft:Fs)';
        mag = 20*log10(abs(fft(fullIR, Nfft)));
        plot(f,mag,'--r','LineWidth',2);

        set(gca,'xscale','log')
        grid on
        xlim([10 Fs/2]);
        % ylim([-80 -30]);
        xlabel('Frequency (Hz)');
        ylabel('Magnitude (dB)');
        
        % plot freefield impulse response
        figure('Name','FF IR','NumberTitle','off','WindowStyle','docked');
        hold on
        box on
        plot(FFIR(:,1),'g')
        plot(FFIR(:,2),'r')
        
        % plot inverese filter for freefield impulse response
        figure('Name','FF inv filter IR','NumberTitle','off','WindowStyle','docked');
        hold on
        box on
        plot(invh(:,1),'g')
        plot(invh(:,2),'r')
    end
    
    %% filter all windowed hrirs
    for i = 1:length(irBank)
        irBank(i).hrirLeft = conv(irBank(i).winIrLeft, invh(:,1));
        irBank(i).hrirRight = conv(irBank(i).winIrRight, invh(:,2));   
    end
    
    %% cut equalized IRs
    figure('Name','raw IR truncation','NumberTitle','off','WindowStyle','docked');
    hold on
    cut_start = 1;
    cut_end = cut_start + 255 + 256;
    for i = 1:length(irBank)
        plot(irBank(i).hrirLeft,'b')
        plot(irBank(i).hrirRight,'r')
    end
    xline(cut_start,'--k')
    xline(cut_end,'--k')
    
    for i = 1:length(irBank)
        irBank(i).hrirLeft = irBank(i).hrirLeft(cut_start:cut_end);
        irBank(i).hrirRight = irBank(i).hrirRight(cut_start:cut_end);
    end
    
end

function IRbank = dfeHRIRs(IRbank, dfe_enabled, plotting, save_fig_folder)
    Fs = unique([IRbank.Fs]);
    
    % remove the ff measurement
    IRbank(isnan([IRbank.azimuth])) = [];
    
    % calculate weights based on solid angle of each cell
    azel = [[IRbank.azimuth]',[IRbank.elevation]'];
    azel = azel + randn(size(azel))*0.001; % add some noise to smooth out plotting
    [xyz(1,:), xyz(2,:), xyz(3,:)] = sph2cart(deg2rad(azel(:,1)),deg2rad(azel(:,2)),1);
    xyz = xyz ./ sqrt(sum(xyz.^2,1));
    [~, ~, voronoiboundary, s] = voronoisphere(xyz);
    
    sn = s./sum(s);
    
    % calculate average magnitude for left & right
    Nfft = 4096;
    mag_ir_avgL = zeros(Nfft,1);
    mag_ir_avgR = zeros(Nfft,1);
    
    for i = 1:length(IRbank)  
        mag_ir_avgL = mag_ir_avgL + (abs(fft(IRbank(i).hrirLeft, Nfft)).^2) * sn(i);
        mag_ir_avgR = mag_ir_avgR + (abs(fft(IRbank(i).hrirRight, Nfft)).^2) * sn(i);
    end
    mag_ir_avgL = sqrt(mag_ir_avgL);
    mag_ir_avgR = sqrt(mag_ir_avgR);

%     % flatten magnitude above 16kHz
%     f = (Fs/Nfft:Fs/Nfft:Fs);
%     mag_ir_avgL = flattenMagHF(f,mag_ir_avgL);
%     mag_ir_avgR = flattenMagHF(f,mag_ir_avgR);
    
    % back to time domain
    ir_avgL = ifft(mag_ir_avgL,'symmetric');
    ir_avgR = ifft(mag_ir_avgR,'symmetric');
    
%     Nfft = 4096/2; % inv fir Nfft and filter length
    Nfft = length(ir_avgL)/2;
    if dfe_enabled   
        % create dfe filters
        dfeL = invFIR('minphase', ir_avgL, Nfft, 0, Nfft, [10 Fs/2], [120 0], 1, Fs);
        dfeR = invFIR('minphase', ir_avgR, Nfft, 0, Nfft, [10 Fs/2], [120 0], 1, Fs);
    else
        % use "flat" filters
        dfeL = [1; zeros(Nfft-1,1)];
        dfeR = [1; zeros(Nfft-1,1)];
    end

    %% filter hrirs with dfe filters
    for i = 1:length(IRbank)
        IRbank(i).dfeHrirLeft = conv(IRbank(i).hrirLeft, dfeL);
        IRbank(i).dfeHrirRight = conv(IRbank(i).hrirRight, dfeR);  
    end
    
    %% cut equalized IRs
    figure('Name','dfe IR truncation','NumberTitle','off','WindowStyle','docked');
    hold on
    cut_start = 1;
    cut_end = cut_start + 255;
    for i = 1:length(IRbank)
        plot(IRbank(i).dfeHrirLeft,'b')
        plot(IRbank(i).dfeHrirRight,'r')
    end
    xline(cut_start,'--k')
    xline(cut_end,'--k')
    
    for i = 1:length(IRbank)
        IRbank(i).dfeHrirLeft = IRbank(i).dfeHrirLeft(cut_start:cut_end);
        IRbank(i).dfeHrirRight = IRbank(i).dfeHrirRight(cut_start:cut_end);
    end
    
%     % calculate target magnitude
%     for i = 1:length(IRbank)
%         for selected_channel = 1:2
%             hrir = IRbank(i).dfehrir(:,selected_channel);
%             Fs = IRbank(i).Fs;
%             Nfft = size(hrir,1)*2;
%             f = (Fs/Nfft:Fs/Nfft:Fs)';
%             mag_target = abs(fft(hrir,Nfft));
%             IRbank(i).f = f(1:Nfft/2);
%             IRbank(i).mag_target(:,selected_channel) = mag_target(1:Nfft/2);
%         end
%     end
    
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
        xlim([100 Fs/2]);
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

function plotMagnitudes(IRbank, type, save_fig_folder)
    %% plot full ir spectrum
    figure('Name',[type ' magnitudes'],'NumberTitle','off','WindowStyle','docked');
    tiledlayout(ceil(sqrt(length(IRbank))),floor(sqrt(length(IRbank))))
    for i = 1:length(IRbank)
        nexttile
        hold on
        Fs = IRbank(i).Fs;
        Nfft = 4096;
        if strcmp(type,'1-measured')
            % left
            fullIR = IRbank(i).fullIrLeft;
            f = (Fs/Nfft:Fs/Nfft:Fs)';
            mag = 20*log10(abs(fft(fullIR, Nfft)));
            plot(f,mag,'-b','LineWidth',1);
            % right
            fullIR = IRbank(i).fullIrRight;
            f = (Fs/Nfft:Fs/Nfft:Fs)';
            mag = 20*log10(abs(fft(fullIR, Nfft)));
            plot(f,mag,'-r','LineWidth',1);
        elseif strcmp(type,'2-win')
            fullIR = IRbank(i).winIrLeft;
            f = (Fs/Nfft:Fs/Nfft:Fs)';
            mag = 20*log10(abs(fft(fullIR, Nfft)));
            plot(f,mag,'-b','LineWidth',1);
            
            fullIR = IRbank(i).winIrRight;
            f = (Fs/Nfft:Fs/Nfft:Fs)';
            mag = 20*log10(abs(fft(fullIR, Nfft)));
            plot(f,mag,'-r','LineWidth',1);
        elseif strcmp(type,'3-raw')
            fullIR = IRbank(i).hrirLeft;
            f = (Fs/Nfft:Fs/Nfft:Fs)';
            mag = 20*log10(abs(fft(fullIR, Nfft)));
            plot(f,mag,'-b','LineWidth',1);

            fullIR = IRbank(i).hrirRight;
            f = (Fs/Nfft:Fs/Nfft:Fs)';
            mag = 20*log10(abs(fft(fullIR, Nfft)));
            plot(f,mag,'-r','LineWidth',1);
        elseif strcmp(type,'4-dfe')
            fullIR = IRbank(i).dfeHrirLeft;
            f = (Fs/Nfft:Fs/Nfft:Fs)';
            mag = 20*log10(abs(fft(fullIR, Nfft)));
            plot(f,mag,'-b','LineWidth',1);

            fullIR = IRbank(i).dfeHrirRight;
            f = (Fs/Nfft:Fs/Nfft:Fs)';
            mag = 20*log10(abs(fft(fullIR, Nfft)));
            plot(f,mag,'-r','LineWidth',1);
        end
        
%         legend('Left Ear','Right Ear','location','southwest')

%         xlim([100 20000]);
        xlim([20 Fs/2]);
        ylim([-40 40]);
        xlabel('Frequency (Hz)');
        ylabel('Magnitude (dB)');
        title(['az: ' num2str(IRbank(i).azimuth) ', el: ' num2str(IRbank(i).elevation) ', dist: ' num2str(IRbank(i).distance)]);
        set(gca,'xscale','log')
    end

    % save figure
    figlen = 16;
    width = 4*figlen;
    height = 3*figlen;
    set(gcf,'PaperPosition',[0 0 width height],'PaperSize',[width height]);
    saveas(gcf,[save_fig_folder 'measured_mag_' type '.jpg'])
%     close

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
        audiowrite([subjectdir '/ambix_wav/' filename], gain*[IRbank(i).dfeHrirLeft IRbank(i).dfeHrirRight], IRbank(i).Fs)
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
       hrirs(i,1,:) = IRbank(i).dfeHrirLeft;
       hrirs(i,2,:) = IRbank(i).dfeHrirRight;
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
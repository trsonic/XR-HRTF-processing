function IRbank = winIRs(IRbank, plotting, save_fig_folder)
    % define window    
    win1 = hann(80).^4;
    win1 = win1(1:end/2);
    win2 = hann(200).^4;
    win2 = win2(end/2+1:end);
    win = [win1; ones(40,1); win2;];
    winshift = length(win1); % how many samples the window should be shifted forward from the peak

    Nwin = length(win);
    if mod(Nwin,2) ~= 0
        disp('wrong window length, must be even')
    end
    
    for i = 1:length(IRbank)
        irLeft = [zeros(Nwin,1); IRbank(i).fullIR(:,1);];
        irRight = [zeros(Nwin,1); IRbank(i).fullIR(:,2);];
        Fs = IRbank(i).Fs;
        
        % get ITD, direct sound sample indices, and direct sound delay
        [IRbank(i).ITD,maxL,maxR,IRbank(i).dlyL,IRbank(i).dlyR] = getITD(irLeft,irRight,Fs);

        % apply window
        winstart = fix(maxL - winshift);
        winend = winstart + Nwin - 1;
        irLeft(1:winstart-1) = 0;
        irLeft(winstart:winend) = irLeft(winstart:winend) .* win;
        irLeft(winend+1:end) = 0;
        
        winstart = fix(maxR - winshift);
        winend = winstart + Nwin - 1;
        irRight(1:winstart-1) = 0;
        irRight(winstart:winend) = irRight(winstart:winend) .* win;
        irRight(winend+1:end) = 0;
        
        % mean time of arrival for both ear signals expressed in samples
        toasmp = round(mean([maxL, maxR]));
        
        % cut irs
        preSamples = fix(0.00075 * Fs);
        afterSamples = 512 - preSamples;
        IRbank(i).winIR(:,1) = irLeft(toasmp-preSamples+1:toasmp+afterSamples);
        IRbank(i).winIR(:,2) = irRight(toasmp-preSamples+1:toasmp+afterSamples);
        
        IRbank(i).toasmp = toasmp - Nwin; % because Nwin of zeroes was added at the beginning
        IRbank(i).maxL = maxL - Nwin;
        IRbank(i).maxR = maxR - Nwin;
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
            IRbank(i).gain = 20*log10(gain_lin);
            IRbank(i).winIR = IRbank(i).winIR * gain_lin;
        end
    end
    
    if strcmp(plotting,'true')
        figure('Name','ETC + win','NumberTitle','off','WindowStyle','docked');
        range = [-1 1];
        mkdir(save_fig_folder)
        for i = 1:length(IRbank) % this will produce lots of graphs
%         for i = randperm(length(IRbank),20) % 20 random directions
            subplot(2,2,1)
            hold on
            yyaxis left
            cla
            plot(IRbank(i).fullIR(:,1),'-b')
            xline(IRbank(i).maxL, '--k','linewidth',2)
%             xline(IRbank(i).toasmp, '--g')
            xlabel('Samples')
            ylabel('Amplitude')
            ylim(range)
            yyaxis right
            cla
            plot(IRbank(i).maxL-winshift:IRbank(i).maxL-winshift+Nwin-1,win,'--')
            ylabel('Linear gain')
            xlim([IRbank(i).maxL-Nwin IRbank(i).maxL+Nwin])
%             title(['Full BRIR (left)' ' azi ' num2str(IRbank(i).azimuth) ' ele ' num2str(IRbank(i).elevation)])
            box on
            
            legend('IR','Peak','Windowing function','Location','NorthWest')

            subplot(2,2,3)
            cla
            hold on
            plot(IRbank(i).winIR(:,1),'-b')
            xlabel('Samples')
            ylabel('Amplitude')
            ylim(range);
            xlim([0 size(IRbank(i).winIR,1)])
%             title(['Windowed BRIR (left)' ' azi ' num2str(IRbank(i).azimuth) ' ele ' num2str(IRbank(i).elevation)])
            box on

            subplot(2,2,2)
            hold on
            yyaxis left
            cla
            plot(IRbank(i).fullIR(:,2),'-b')
            xline(IRbank(i).maxR, '--k','linewidth',2)
            xline(IRbank(i).toasmp, '--g')
            xlabel('Samples')
            ylabel('Amplitude')
            ylim(range)
            yyaxis right
            cla
            plot(IRbank(i).maxR-winshift:IRbank(i).maxR-winshift+Nwin-1,win,'--')
            ylabel('Linear gain')
            xlim([IRbank(i).maxR-Nwin IRbank(i).maxR+Nwin])
%             title('Right Raw HRIR')
            box on


            subplot(2,2,4)
            cla
            hold on
            plot(IRbank(i).winIR(:,2),'-b')
            xlabel('Samples')
            ylabel('Amplitude')
            ylim(range);
            xlim([0 size(IRbank(i).winIR,1)])
%             title('Right Windowed HRIR')
            box on

            
            % save figure
            figlen = 8;
            width = 4*figlen;
            height = 2*figlen;
            set(gcf,'Units','centimeters','PaperPosition',[0 0 width height],'PaperSize',[width height]);
%             saveas(gcf,[save_fig_folder 'windowing_fig' num2str(i)  '.jpg'])
            saveas(gcf,[save_fig_folder 'windowing_fig' num2str(i)  '.pdf'])
            close
        end
    end
end
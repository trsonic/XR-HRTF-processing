function irBank = normalizeIRs(irBank, plotting, save_fig_folder)
    % find the free-field measurements and calculate inverse filters
    for i = find([irBank.ref])
        Fs = irBank(i).Fs;
        irBank(i).invh(:,1) = createInverseFilter(irBank(i).winIR(:,1), Fs, 12, [0 0]);
        irBank(i).invh(:,2) = createInverseFilter(irBank(i).winIR(:,2), Fs, 12, [0 0]);
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
%             grid on
            xlim([20 Fs/2]);
            ylim([-35 35]);
            legend('Left channel', 'Right channel', 'Left inverse filter', 'Right inverse filter','location','NorthEast')
            xlabel('Frequency (Hz)');
            ylabel('Magnitude (dB)');
            
            % save figure
            figlen = 4;
            width = 4*figlen;
            height = 2.5*figlen;
            set(gcf,'Units','centimeters','PaperPosition',[0 0 width height],'PaperSize',[width height]);
            mkdir(save_fig_folder)
%             saveas(gcf,[save_fig_folder 'normalization-' num2str(i) '.png'])
            saveas(gcf,[save_fig_folder 'normalization-' num2str(i) '.pdf']) 
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
        plotting = 'false';
%         if irBank(i).azimuth == 90 && irBank(i).elevation == 0
%             plotting = 'true';
%         end
        dist = 1.5; % reference distance
        hd = 0.16; % head diameter / ear to ear distance
        dd = dist + (hd/2) * sin(deg2rad(-irBank(i).azimuth)) * cos(deg2rad(irBank(i).elevation));
        lfe_amp = 20*log10(dd/dist);
        irBank(i).rawHRIR(:,1) = LFextension(irBank(i).rawHRIR(:,1), Fs, lfe_amp, plotting);
        irBank(i).rawHRIR(:,2) = LFextension(irBank(i).rawHRIR(:,2), Fs, -lfe_amp, plotting);
    end
    
    %% cut equalized IRs
    figure('Name','raw IR truncation','NumberTitle','off','WindowStyle','docked');
    hold on
    cut_start = 1;
    cut_end = 300;
    for i = 1:length(irBank)
        plot(irBank(i).rawHRIR(:,1),'b')
        plot(irBank(i).rawHRIR(:,2),'r')
    end
    xline(cut_start,'--k')
    xline(cut_end,'--k')
    
    for i = 1:length(irBank)
        irBank(i).rawHRIR = irBank(i).rawHRIR(cut_start:cut_end,:);
    end


function lfe_h = LFextension(h, Fs, lfe_amp, plotting)    
    [acor, lag] = xcorr(h,minph(h));
    [~,index] = max(acor);
    shift = lag(index);

    lfe_kronecker = zeros(length(h),1);
    lfe_kronecker(shift) = 1 * 10^(lfe_amp/20);
    
    xfreq = 250; % crossover frequency
    filter_order = 2;    
    [B_highpass, A_highpass] = butter( filter_order, xfreq/Fs*2, 'high' );            
    [B_lowpass,  A_lowpass ] = butter( filter_order, xfreq/Fs*2, 'low'  );
    output_low = filter(B_lowpass, A_lowpass, lfe_kronecker);
    output_low = filter(B_lowpass, A_lowpass, output_low);
    output_high = filter(B_highpass, A_highpass, h);
    output_high = filter(B_highpass, A_highpass, output_high);
        
    lfe_h = output_low + output_high;
    
    if strcmp(plotting,'true')
        %% simple plot
        figure('Name','lf extension','NumberTitle','off','WindowStyle','docked');
        hold on
        [f,mag] = getMagnitude(h,Fs,'log');
        plot(f,mag,'-b','LineWidth',1);
        [f,mag] = getMagnitude(output_low,Fs,'log');
        plot(f,mag,'--g','LineWidth',2);
        [f,mag] = getMagnitude(output_high,Fs,'log');
        plot(f,mag,'--r','LineWidth',2);
        [f,mag] = getMagnitude(lfe_h,Fs,'log');
        plot(f,mag,'-m','LineWidth',1);
        legend('Original','Low-passed','High-passed','Extended','location','northwest')
        xlim([20 Fs/2]);
        ylim([-30 20]);
        set(gca,'xscale','log')
        xlabel('Frequency (Hz)');
        ylabel('Magnitude (dB)');
        box on
        
        % save figure
        figlen = 4;
        width = 4*figlen;
        height = 2.5*figlen;
        set(gcf,'Units','centimeters','PaperPosition',[0 0 width height],'PaperSize',[width height]);
        saveas(gcf,'lfe.pdf')  
        
%         %% advanced plot
%         figure('Name','lf extension','NumberTitle','off','WindowStyle','docked');
%         subplot(2,2,1)
%         hold on
%         plot(h)
%         plot(lfe_kronecker)
%         plot(output_low,'--','LineWidth',1)
%         plot(output_high,'--','LineWidth',1)
%         plot(lfe_h,'-','LineWidth',2)
%         xline(shift,'--k')
% 
%         legend('original','kronecker','low','high','extended','location','northeast')
%         xlim([0 256]);
%         ylim([-1 1]);
% 
%         subplot(2,2,2)
%         hold on
%         [gd,f] = grpdelay(h,1,2^10,'whole',Fs);
%         plot(f,gd,'-b','LineWidth',1);
%         [gd,f] = grpdelay(lfe_h,1,2^10,'whole',Fs);
%         plot(f,gd,'-r','LineWidth',1);
%         xlim([20 Fs/2]), ylim([-200 500])
%         xlabel('Frequency (Hz)'), ylabel('group delay in samples')
%         set(gca,'xscale','log')
%         legend('original','extended','location','southwest')
% 
%         subplot(2,2,3)
%         hold on
%         [f,mag] = getMagnitude(h,Fs,'log');
%         plot(f,mag,'-b','LineWidth',1);
%         [f,mag] = getMagnitude(output_low,Fs,'log');
%         plot(f,mag,'--g','LineWidth',2);
%         [f,mag] = getMagnitude(output_high,Fs,'log');
%         plot(f,mag,'--r','LineWidth',2);
%         [f,mag] = getMagnitude(lfe_h,Fs,'log');
%         plot(f,mag,'-m','LineWidth',1);
%         legend('original','low','high','extended','location','northwest')
%         xlim([20 Fs/2]);
%         ylim([-30 20]);
%         set(gca,'xscale','log')
%         xlabel('Frequency (Hz)');
%         ylabel('Magnitude (dB)');
% 
%         subplot(2,2,4)
%         hold on
%         [f,mag] = getPhase(h,Fs);
%         plot(f,mag,'-b','LineWidth',1);
%         [f,mag] = getPhase(output_low,Fs);
%         plot(f,mag,'--g','LineWidth',2);
%         [f,mag] = getPhase(output_high,Fs);
%         plot(f,mag,'--r','LineWidth',2);
%         [f,mag] = getPhase(lfe_h,Fs);
%         plot(f,mag,'-m','LineWidth',1);
%         legend('original','low','high','extended','location','southwest')
%         xlim([20 Fs/2]);
%     %     ylim([-40 40]);
%         set(gca,'xscale','log')
%         xlabel('Frequency (Hz)');
%         ylabel('Phase');
    end
end
end
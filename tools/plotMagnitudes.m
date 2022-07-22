function plotMagnitudes(irBank, type, save_fig_folder)
    plots_num = 51;
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

%             legend('Left Ear','Right Ear','location','southwest')
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
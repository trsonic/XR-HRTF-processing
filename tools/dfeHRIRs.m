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
    for i = 1:length(irBank)
        if i == 1
            Nfft = size(irBank(i).rawHRIR,1);
            mag_ir_avgL = zeros(Nfft,1);
            mag_ir_avgR = zeros(Nfft,1);
        end
        mag_ir_avgL = mag_ir_avgL + (abs(fft(irBank(i).rawHRIR(:,1), Nfft)).^2) * sn(i);
        mag_ir_avgR = mag_ir_avgR + (abs(fft(irBank(i).rawHRIR(:,2), Nfft)).^2) * sn(i);
    end

    mag_ir_avgLR = sqrt(mag_ir_avgL/2 + mag_ir_avgR/2);
    mag_ir_avgL = sqrt(mag_ir_avgL);
    mag_ir_avgR = sqrt(mag_ir_avgR);

    % back to time domain
    ir_avgLR = ifft(mag_ir_avgLR,'symmetric');
    ir_avgLR = circshift(ir_avgLR,Nfft/2,1);

    ir_avgL = ifft(mag_ir_avgL,'symmetric');
    ir_avgR = ifft(mag_ir_avgR,'symmetric');
    ir_avgL = circshift(ir_avgL,Nfft/2,1);
    ir_avgR = circshift(ir_avgR,Nfft/2,1);
    
    % inv fir Nfft and filter length
    Nfft = length(ir_avgL);
    if dfe_enabled   
        % create dfe filters
        dfeLR = createInverseFilter(ir_avgLR, Fs, 12, [0 0]);
    else
        % use "flat" filters
        dfeLR = [1; zeros(Nfft-1,1)];
    end

    %% filter hrirs with dfe filters
    for i = 1:length(irBank)
        irBank(i).dfeHRIR(:,1) = conv(dfeLR,irBank(i).rawHRIR(:,1));
        irBank(i).dfeHRIR(:,2) = conv(dfeLR,irBank(i).rawHRIR(:,2));  
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
%         grid on
        [f,mag] = getMagnitude(ir_avgLR,Fs,'log');
        plot(f,mag,'-b','LineWidth',1);
        [f,mag] = getMagnitude(ir_avgL,Fs,'log');
        plot(f,mag,'-g','LineWidth',1);
        [f,mag] = getMagnitude(ir_avgR,Fs,'log');
        plot(f,mag,'-r','LineWidth',1);
        [f,mag] = getMagnitude(dfeLR,Fs,'log');
        plot(f,mag,'--b','LineWidth',1.5);
                
        set(gca,'xscale','log')
        xlim([100 Fs/2]);
        ylim([-20 20]);
        ylabel('Magnitude (dB)')
        xlabel('Frequency (Hz)')
        legend('Left and Right Ear Average', 'Left Ear Average', 'Right Ear Average', 'DFE Filter','Location','southwest')
%         title('DFE filter')
        
        % save figure
        figlen = 4;
        width = 4*figlen;
        height = 2.5*figlen;
        set(gcf,'Units','centimeters','PaperPosition',[0 0 width height],'PaperSize',[width height]);
        mkdir(save_fig_folder)
        saveas(gcf,[save_fig_folder 'dfe.pdf'])        
%         close
    end
end
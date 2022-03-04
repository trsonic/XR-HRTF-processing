close all
clear

subjectdir = 'data/20211126-XR-TR/';
% subjectdir = 'data/20211126-XR-Gavin/';
load([subjectdir 'irBankProcessed.mat'])

% does barycentric interpolation of hrirs
addpath('../API_MO/API_MO/')

addpath('tools/TriangleRayIntersection/')
addpath('tools/')

interpHrirBank = interpHRIRs(irBank,'raw');
saveAsSofa(interpHrirBank,subjectdir,'raw')

interpHrirBank = interpHRIRs(irBank,'dfe');
saveAsSofa(interpHrirBank,subjectdir,'dfe')

function interpHrirBank = interpHRIRs(hrirBank, type)
    % get minimum phase hrirs
    for i = 1:length(hrirBank)
        if strcmp(type,'raw')
            left_hrir(i,:) = minph(hrirBank(i).rawHRIR(:,1));
            right_hrir(i,:) = minph(hrirBank(i).rawHRIR(:,2));
        elseif strcmp(type,'dfe')
            left_hrir(i,:) = minph(hrirBank(i).dfeHRIR(:,1));
            right_hrir(i,:) = minph(hrirBank(i).dfeHRIR(:,2));
        end

        ITD(i) = hrirBank(i).ITD;
    end

    % get interpolated dirs matrix
    azel = [[hrirBank.azimuth]' [hrirBank.elevation]'];
    azel_interp = [];
%     for az = -180:2:180
%         for el = -90:2:90
%             azel_interp = [azel_interp; az el];
%         end
%     end
    ls = getLebedevSphere(2354);
    dirs = [];
    [dirs(:,1), dirs(:,2), ~] = cart2sph(ls.x,ls.y,ls.z);
    azel_interp = rad2deg(dirs);

    % get barycentric weights
    [bid, bw] = barycentric_interpolation(azel, azel_interp);

    % calculate interpolated minimum phase hrirs
    for i = 1:length(azel_interp)
        interpHrirBank(i).azimuth = azel_interp(i,1);
        interpHrirBank(i).elevation = azel_interp(i,2);
        interpHrirBank(i).Fs = unique([hrirBank.Fs]);
        interpHrirBank(i).ITD =       (ITD(bid(i,1))*bw(i,1) + ...
                                       ITD(bid(i,2))*bw(i,2) + ...
                                       ITD(bid(i,3))*bw(i,3)) / ...
                                       sum(bw(i,:));
        interpHrirBank(i).left_hrir = (left_hrir(bid(i,1),:)*bw(i,1) + ...
                                       left_hrir(bid(i,2),:)*bw(i,2) + ...
                                       left_hrir(bid(i,3),:)*bw(i,3)) / ...
                                       sum(bw(i,:));

        interpHrirBank(i).right_hrir = (right_hrir(bid(i,1),:)*bw(i,1) + ...
                                       right_hrir(bid(i,2),:)*bw(i,2) + ...
                                       right_hrir(bid(i,3),:)*bw(i,3)) / ...
                                       sum(bw(i,:));
        

        
        if strcmp(type,'raw')
            interpHrirBank(i).rawHRIR = [interpHrirBank(i).left_hrir; interpHrirBank(i).right_hrir]';                        
            interpHrirBank(i).rawHRIR = injectITD(interpHrirBank(i).rawHRIR,interpHrirBank(i).ITD,interpHrirBank(i).Fs);
        elseif strcmp(type,'dfe')
            interpHrirBank(i).dfeHRIR = [interpHrirBank(i).left_hrir; interpHrirBank(i).right_hrir]';                        
            interpHrirBank(i).dfeHRIR = injectITD(interpHrirBank(i).dfeHRIR,interpHrirBank(i).ITD,interpHrirBank(i).Fs);
        end

    end

    plotHMFmags(interpHrirBank)
    plotITD([interpHrirBank.azimuth],[interpHrirBank.elevation],[interpHrirBank.ITD], [-1000 1000])
end

function [idx, bweights] = barycentric_interpolation(azel, azel_interp)
    [vx(:,1), vx(:,2), vx(:,3)] = sph2cart(deg2rad(azel(:,1)),deg2rad(azel(:,2)),1);
    [vqx(:,1), vqx(:,2), vqx(:,3)] = sph2cart(deg2rad(azel_interp(:,1)),deg2rad(azel_interp(:,2)),1);

    faces = convhull(vx, 'simplify', true);
    vert1 = vx(faces(:,1),:);
    vert2 = vx(faces(:,2),:);
    vert3 = vx(faces(:,3),:);
    
    for i = 1:size(vqx,1)
        epsilon = 1e-5;
        [intersect, ~, bary_weights_1, bary_weights_2, ~] = TriangleRayIntersection([0 0 0], vqx(i,:), vert1, vert2, vert3,'lineType','ray','planeType','two sided','border','inclusive','eps',epsilon);
        intersect = find(intersect);       % indices of the hit triangles
        assert( ~isempty(intersect), 'intersect is empty');  % check if any hit points have been found

        % pick the first hit triangle
        first_face_index = intersect(1);
        id = faces(first_face_index, :);
        w1 = bary_weights_1(first_face_index);
        w2 = bary_weights_2(first_face_index);
        w3 = 1 - w1 - w2;
        bw = [w3 w1 w2];
        assert( sum(bw) <= 1+epsilon && sum(bw) >= 1-epsilon && bw(1) >= -epsilon && bw(2) >= -epsilon && bw(3) >= -epsilon, 'barycentric coordinates error');
        
        idx(i,:) = id;
        bweights(i,:) = bw;
    end
end

function h_ITD = injectITD(h,ITD,Fs)
    shift_l = (750 - ITD/2) * 10^-6 * Fs;
    shift_r = (750 + ITD/2) * 10^-6 * Fs;
    h_l = fraccircshift(h(:,1),shift_l);
    h_r = fraccircshift(h(:,2),shift_r);
    h_ITD = [h_l h_r];
    
%     % plot
%     subplot(2,1,1)
%     hold on
%     plot(h)
%     subplot(2,1,2)
%     hold on
%     plot(h_ITD)
end

function plotHMFmags(irBank)
    zlimit = [-40 20];
    %% plot surface plots
    figure('Name','Surf plots','NumberTitle','off','WindowStyle','docked')
    tiledlayout(3,2)

    % horizontal plane
    az = -180:1:180;
    el = 0;
    left_hrirs = [];
    right_hrirs = [];
    for i = 1:length(az)
        dist = distance([irBank.elevation], [irBank.azimuth], el, az(i));
        [~, idx] = min(dist);
        left_hrirs(:,i) = irBank(idx).left_hrir;
        right_hrirs(:,i) = irBank(idx).right_hrir;
    end
    nexttile
    magsurf(irBank(idx).Fs,az,left_hrirs,'Azimuth (deg)','Horizontal plane - left ear',zlimit)
    nexttile
    magsurf(irBank(idx).Fs,az,right_hrirs,'Azimuth (deg)','Horizontal plane - right ear',zlimit)

    % median plane
    az = 0;
    el = -180:1:180;
    left_hrirs = [];
    right_hrirs = [];
    for i = 1:length(el)
        dist = distance([irBank.elevation], [irBank.azimuth], el(i), az);
        [~, idx] = min(dist);
        left_hrirs(:,i) = irBank(idx).left_hrir;
        right_hrirs(:,i) = irBank(idx).right_hrir;
    end
    nexttile
    magsurf(irBank(idx).Fs,el,left_hrirs,'Elevation (deg)','Median plane - left ear',zlimit)
    nexttile
    magsurf(irBank(idx).Fs,el,right_hrirs,'Elevation (deg)','Median plane - right ear',zlimit)

    % frontal plane
    az = 90;
    el = -180:1:180;
    left_hrirs = [];
    right_hrirs = [];
    for i = 1:length(el)
        dist = distance([irBank.elevation], [irBank.azimuth], el(i), az);
        [~, idx] = min(dist);
        left_hrirs(:,i) = irBank(idx).left_hrir;
        right_hrirs(:,i) = irBank(idx).right_hrir;
    end
    nexttile
    magsurf(irBank(idx).Fs,el,left_hrirs,'Elevation (deg)','Frontal plane - left ear',zlimit)
    nexttile
    magsurf(irBank(idx).Fs,el,right_hrirs,'Elevation (deg)','Frontal plane - right ear',zlimit)
end

function magsurf(Fs,ang,hrir,ytit,tit,zlimit)
    Nfft = 1024;
    f = ((0:Nfft-1)*Fs/Nfft)';
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

function plotITD(az, el, val, zlimit)
    figure('Name','ITD','NumberTitle','off','WindowStyle','docked');
    azimuth = -180:1:180;
    elevation = -90:1:90;

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
    zlim(zlimit)
    caxis(zlimit)
    xlabel('Azimuth (deg)')
    ylabel('Elevation (deg)')
%     title('MFSD')
    colorbar
end

function saveAsSofa(IRbank, subjectdir, type)
    % Start SOFA
    SOFAstart
    
    % Get an empy conventions structure
    Obj = SOFAgetConventions('SimpleFreeFieldHRIR');
    
    % remove the ff measurement
    IRbank(isnan([IRbank.azimuth])) = [];

    for i = 1:length(IRbank)
        if strcmp(type,'raw')
            hrirs(i,:,:) = IRbank(i).rawHRIR';
        elseif strcmp(type,'dfe')
            hrirs(i,:,:) = IRbank(i).dfeHRIR';
        end
       azi = IRbank(i).azimuth;
       ele = IRbank(i).elevation;
       dist = 1.5;

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
    SOFAfn=fullfile(subjectdir,['xr-hrtf-interp-' type '.sofa']);
    disp(['Saving:  ' SOFAfn]);
    SOFAsave(SOFAfn, Obj, compression);
end
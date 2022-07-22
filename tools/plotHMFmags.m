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
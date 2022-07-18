function plotAzEl(az, el, val, lim)
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
    zlim(lim)
    caxis(lim)
    xlabel('Azimuth (deg)')
    ylabel('Elevation (deg)')
%     title('MFSD')
    colorbar
end
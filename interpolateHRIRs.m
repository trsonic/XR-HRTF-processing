close all
clear

%% add paths and initialize SOFA API
addpath('../API_MO/API_MO/')
addpath('tools/TriangleRayIntersection/')
SOFAstart;

sofafile = 'data/20201217-122pt-2.5m-dayton_vt/xr_122pt.sofa';

%% load sofa file
hrtf = SOFAload(sofafile);
hrirBank = [];
apparentSourceVector = SOFAcalculateAPV(hrtf); % Calculate the source position from a listener point of view
for i = 1:length(apparentSourceVector)
    hrirBank(i).azimuth = apparentSourceVector(i,1);
    hrirBank(i).elevation = apparentSourceVector(i,2);
    hrirBank(i).left_hrir = squeeze(hrtf.Data.IR(i,1,:))';
    hrirBank(i).right_hrir = squeeze(hrtf.Data.IR(i,2,:))';
    hrirBank(i).Fs = hrtf.Data.SamplingRate;        
end

% get magnitudes
for i = 1:length(hrirBank)
    Fs = hrirBank(i).Fs;
    Nfft = size(hrirBank(i).left_hrir,2)*2;
    f = (Fs/Nfft:Fs/Nfft:Fs);
    hrirBank(i).left_mag = abs(fft(hrirBank(i).left_hrir,Nfft));
    hrirBank(i).right_mag = abs(fft(hrirBank(i).right_hrir,Nfft));
    hrirBank(i).left_mag = hrirBank(i).left_mag(1:Nfft/2);
    hrirBank(i).right_mag = hrirBank(i).right_mag(1:Nfft/2);
    hrirBank(i).f = f(1:Nfft/2);
    
    left_mag(i,:) = hrirBank(i).left_mag;
    right_mag(i,:) = hrirBank(i).right_mag;
end

% get interpolated dirs matrix
azel = [[hrirBank.azimuth]' [hrirBank.elevation]'];
azel_interp = [];
for az = -180:2:180
    for el = -90:2:90
        azel_interp = [azel_interp; az el];
    end
end

% get barycentric weights
[bid, bw] = barycentric_interpolation(azel, azel_interp);

% calculate interpolated magnitudes
for i = 1:length(azel_interp)
    interpHrirBank(i).azimuth = azel_interp(i,1);
    interpHrirBank(i).elevation = azel_interp(i,2);
    interpHrirBank(i).f = unique([hrirBank.f]);
    interpHrirBank(i).left_mag = (left_mag(bid(i,1),:)*bw(i,1) + ...
                                   left_mag(bid(i,2),:)*bw(i,2) + ...
                                   left_mag(bid(i,3),:)*bw(i,3)) / ...
                                   sum(bw(i,:));
                               
    interpHrirBank(i).right_mag = (right_mag(bid(i,1),:)*bw(i,1) + ...
                                   right_mag(bid(i,2),:)*bw(i,2) + ...
                                   right_mag(bid(i,3),:)*bw(i,3)) / ...
                                   sum(bw(i,:));
end

%% plot surface plots
figure('Name','Surf plots','NumberTitle','off','WindowStyle','docked')
tiledlayout(3,2)

% horizontal plane
az = -180:1:180;
el = 0;
left_mags = [];
right_mags = [];
for i = 1:length(az)
    dist = distance([interpHrirBank.elevation], [interpHrirBank.azimuth], el, az(i));
    [~, idx] = min(dist);
    left_mags(:,i) = interpHrirBank(idx).left_mag;
    right_mags(:,i) = interpHrirBank(idx).right_mag;
end
nexttile
magsurf(az,unique([interpHrirBank.f]),left_mags,'Azimuth (deg)','Horizontal plane - left ear',[-60 20])
nexttile
magsurf(az,unique([interpHrirBank.f]),right_mags,'Azimuth (deg)','Horizontal plane - right ear',[-60 20])

% median plane
az = 0;
el = -180:1:180;
left_mags = [];
right_mags = [];
for i = 1:length(el)
    dist = distance([interpHrirBank.elevation], [interpHrirBank.azimuth], el(i), az);
    [~, idx] = min(dist);
    left_mags(:,i) = interpHrirBank(idx).left_mag;
    right_mags(:,i) = interpHrirBank(idx).right_mag;
end
nexttile
magsurf(el,unique([interpHrirBank.f]),left_mags,'Elevation (deg)','Median plane - left ear',[-60 20])
nexttile
magsurf(el,unique([interpHrirBank.f]),right_mags,'Elevation (deg)','Median plane - right ear',[-60 20])

% frontal plane
az = 90;
el = -180:1:180;
left_mags = [];
right_mags = [];
for i = 1:length(el)
    dist = distance([interpHrirBank.elevation], [interpHrirBank.azimuth], el(i), az);
    [~, idx] = min(dist);
    left_mags(:,i) = interpHrirBank(idx).left_mag;
    right_mags(:,i) = interpHrirBank(idx).right_mag;
end
nexttile
magsurf(el,unique([interpHrirBank.f]),left_mags,'Elevation (deg)','Frontal plane - left ear',[-60 20])
nexttile
magsurf(el,unique([interpHrirBank.f]),right_mags,'Elevation (deg)','Frontal plane - right ear',[-60 20])

%%
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

function magsurf(ang,f,magnitude,ytit,tit,zlimit)
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

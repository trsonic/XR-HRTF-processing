function interpHrirBank = interpHRIRs(hrirBank, type, method)

    %% get hrirs
    for i = 1:length(hrirBank)
        if strcmp(method,'minph')
            if strcmp(type,'raw')
                left_hrir(i,:) = minph(hrirBank(i).rawHRIR(:,1));
                right_hrir(i,:) = minph(hrirBank(i).rawHRIR(:,2));
            elseif strcmp(type,'dfe')
                left_hrir(i,:) = minph(hrirBank(i).dfeHRIR(:,1));
                right_hrir(i,:) = minph(hrirBank(i).dfeHRIR(:,2));
            end
        elseif strcmp(method,'align')
            if strcmp(type,'raw')
                left_hrir(i,:) = hrirBank(i).rawHRIR(:,1);
                right_hrir(i,:) = hrirBank(i).rawHRIR(:,2);
            elseif strcmp(type,'dfe')
                left_hrir(i,:) = hrirBank(i).dfeHRIR(:,1);
                right_hrir(i,:) = hrirBank(i).dfeHRIR(:,2);
            end
        end

        ITD(i) = hrirBank(i).ITD;
    end

    figure('Name','ITD','NumberTitle','off','WindowStyle','docked');
    plotAzElM([hrirBank.azimuth],[hrirBank.elevation],[hrirBank.ITD], [-800 800],'ITD (μs)','contralateral','ipsilateral')

    %% get interpolated dirs matrix
    azel = [[hrirBank.azimuth]' [hrirBank.elevation]'];
    azel_interp = [];
    ls = getLebedevSphere(4334);
    dirs = [];
    [dirs(:,1), dirs(:,2), ~] = cart2sph(ls.x,ls.y,ls.z);
    azel_interp = rad2deg(dirs);

%     %% add HMF points
%     for az = -180:1:180
%         el = 0;
%         azel_interp = [azel_interp; az el];
%     end
%     for el = -90:1:90
%         az = 0;
%         azel_interp = [azel_interp; az el];
%         az = 180;
%         azel_interp = [azel_interp; az el];
%     end
%     for el = -90:1:90
%         az = 90;
%         azel_interp = [azel_interp; az el];
%         az = -90;
%         azel_interp = [azel_interp; az el];
%     end

    %% get barycentric weights
    [bid, bw] = barycentric_interpolation(azel, azel_interp);

    %% calculate interpolated hrirs
    if strcmp(method,'minph')
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

    elseif strcmp(method,'align')
        for i = 1:length(azel_interp)
            interpHrirBank(i).azimuth = azel_interp(i,1);
            interpHrirBank(i).elevation = azel_interp(i,2);
            interpHrirBank(i).Fs = unique([hrirBank.Fs]);
            interpHrirBank(i).ITD =       (ITD(bid(i,1))*bw(i,1) + ...
                                           ITD(bid(i,2))*bw(i,2) + ...
                                           ITD(bid(i,3))*bw(i,3)) / ...
                                           sum(bw(i,:));
    
            s1 = (ITD(bid(i,1)) - interpHrirBank(i).ITD)/2 * 10^-6 * interpHrirBank(i).Fs;
            s2 = (ITD(bid(i,2)) - interpHrirBank(i).ITD)/2 * 10^-6 * interpHrirBank(i).Fs;
            s3 = (ITD(bid(i,3)) - interpHrirBank(i).ITD)/2 * 10^-6 * interpHrirBank(i).Fs;
            interpHrirBank(i).left_hrir = (fraccircshift(left_hrir(bid(i,1),:),s1)*bw(i,1) + ...
                                           fraccircshift(left_hrir(bid(i,2),:),s2)*bw(i,2) + ...
                                           fraccircshift(left_hrir(bid(i,3),:),s3)*bw(i,3)) / ...
                                           sum(bw(i,:));
    
    %         figure('Name','IR Alignment','NumberTitle','off','WindowStyle','docked')
    %         subplot(2,2,1)
    %         hold on
    %         plot(left_hrir(bid(i,1),:))
    %         plot(left_hrir(bid(i,2),:))
    %         plot(left_hrir(bid(i,3),:))
    %         subplot(2,2,3)
    %         hold on
    %         plot(fraccircshift(left_hrir(bid(i,1),:),s1))
    %         plot(fraccircshift(left_hrir(bid(i,2),:),s2))
    %         plot(fraccircshift(left_hrir(bid(i,3),:),s3))
    
            s1 = -(ITD(bid(i,1)) - interpHrirBank(i).ITD)/2 * 10^-6 * interpHrirBank(i).Fs;
            s2 = -(ITD(bid(i,2)) - interpHrirBank(i).ITD)/2 * 10^-6 * interpHrirBank(i).Fs;
            s3 = -(ITD(bid(i,3)) - interpHrirBank(i).ITD)/2 * 10^-6 * interpHrirBank(i).Fs;
            interpHrirBank(i).right_hrir = (fraccircshift(right_hrir(bid(i,1),:),s1)*bw(i,1) + ...
                                           fraccircshift(right_hrir(bid(i,2),:),s2)*bw(i,2) + ...
                                           fraccircshift(right_hrir(bid(i,3),:),s3)*bw(i,3)) / ...
                                           sum(bw(i,:));
    %         subplot(2,2,2)
    %         hold on
    %         plot(right_hrir(bid(i,1),:))
    %         plot(right_hrir(bid(i,2),:))
    %         plot(right_hrir(bid(i,3),:))
    %         subplot(2,2,4)
    %         hold on
    %         plot(fraccircshift(right_hrir(bid(i,1),:),s1))
    %         plot(fraccircshift(right_hrir(bid(i,2),:),s2))
    %         plot(fraccircshift(right_hrir(bid(i,3),:),s3))
            
            if strcmp(type,'raw')
                interpHrirBank(i).rawHRIR = [interpHrirBank(i).left_hrir; interpHrirBank(i).right_hrir]';                        
            elseif strcmp(type,'dfe')
                interpHrirBank(i).dfeHRIR = [interpHrirBank(i).left_hrir; interpHrirBank(i).right_hrir]';                        
            end
        end
    end

    plotHMFmags(interpHrirBank)
%     plotITD([interpHrirBank.azimuth],[interpHrirBank.elevation],[interpHrirBank.ITD], [-1000 1000])
    figure('Name','ITD interpolated','NumberTitle','off','WindowStyle','docked');
    plotAzElM([interpHrirBank.azimuth],[interpHrirBank.elevation],[interpHrirBank.ITD], [-800 800],'ITD (μs)','contralateral','ipsilateral')


    function [idx, bweights] = barycentric_interpolation(azel, azel_interp)
        [vx(:,1), vx(:,2), vx(:,3)] = sph2cart(deg2rad(azel(:,1)),deg2rad(azel(:,2)),1);
        [vqx(:,1), vqx(:,2), vqx(:,3)] = sph2cart(deg2rad(azel_interp(:,1)),deg2rad(azel_interp(:,2)),1);
    
        faces = convhull(vx, 'simplify', true);
        vert1 = vx(faces(:,1),:);
        vert2 = vx(faces(:,2),:);
        vert3 = vx(faces(:,3),:);
        
        for vi = 1:size(vqx,1)
            epsilon = 1e-5;
            [intersect, ~, bary_weights_1, bary_weights_2, ~] = TriangleRayIntersection([0 0 0], vqx(vi,:), vert1, vert2, vert3,'lineType','ray','planeType','two sided','border','inclusive','eps',epsilon);
            intersect = find(intersect);       % indices of the hit triangles
            assert( ~isempty(intersect), 'intersect is empty');  % check if any hit points have been found
    
            % pick the first hit triangle
            first_face_index = intersect(1);
            id = faces(first_face_index, :);
            w1 = bary_weights_1(first_face_index);
            w2 = bary_weights_2(first_face_index);
            w3 = 1 - w1 - w2;
            bwx = [w3 w1 w2];
            assert( sum(bwx) <= 1+epsilon && sum(bwx) >= 1-epsilon && bwx(1) >= -epsilon && bwx(2) >= -epsilon && bwx(3) >= -epsilon, 'barycentric coordinates error');
            
            idx(vi,:) = id;
            bweights(vi,:) = bwx;
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
end
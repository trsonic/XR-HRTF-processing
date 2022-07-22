function saveAsAmbix(IRbank, subjectdir)

    % remove the ff measurement
    IRbank(isnan([IRbank.azimuth])) = [];
    
    % save wav files
    mkdir(subjectdir,'ambix_wav_raw')
    mkdir(subjectdir,'ambix_wav_dfe')
    mkdir(subjectdir,'ambix_wav_hpeq')
    for i = 1:length(IRbank)
        filename = ['azi_' num2str(IRbank(i).azimuth,'%.2f') '_ele_' num2str(IRbank(i).elevation,'%.2f') '.wav'];
        disp(filename);
        IRbank(i).filename = filename;
        gain = 0.5; %10^(-3/20);
        audiowrite([subjectdir '/ambix_wav_raw/' filename], gain * IRbank(i).rawHRIR, IRbank(i).Fs)
        audiowrite([subjectdir '/ambix_wav_dfe/' filename], gain * IRbank(i).dfeHRIR, IRbank(i).Fs)
        rawHRIR = IRbank(i).rawHRIR;
        
        if exist([subjectdir '/hpeq/' 'hpeq.wav']) ~= 0
            [hpeq, Fs] = audioread([subjectdir '/hpeq/' 'hpeq.wav']);
            hpeqHRIR = [conv(rawHRIR(:,1),hpeq(:,1)) conv(rawHRIR(:,2),hpeq(:,2))];
            hpeqHRIR = hpeqHRIR(1:size(rawHRIR,1),:);
            audiowrite([subjectdir '/ambix_wav_hpeq/' filename], gain * hpeqHRIR, IRbank(i).Fs)
        end
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
%         mkdir(subjectdir,'ambix_configs')
%         out_path = [subjectdir 'ambix_configs/xr-hrtf-' num2str(order) 'OA-' layout];
        out_path = [subjectdir 'ambix_wav_dfe/xr-hrtf-' num2str(order) 'OA-' layout];  
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

    function [ val ] = SPKR_ARRAY(azis,eles,dist)
        val.name = 'xr-hrtf';
        % X Y Z
        [Sp(:,1),Sp(:,2),Sp(:,3)] = sph2cart(deg2rad(azis),deg2rad(eles),dist);    
        % spherical coordinate representation
        [val.az, val.el, val.r] = cart2sph(Sp(:,1),Sp(:,2),Sp(:,3));
        
        % unit vectors
        [val.x, val.y, val.z] = sph2cart(val.az, val.el, 1);
       
        for spk = 1:length(azis)
            val.id(spk) = {['spk' num2str(spk)]};
        end
    end
end
%%%%
%%
%% This is the main script for the XR-based HRTF Measurement processing.
%% Author: Tomasz Rudzki (tr837@york.ac.uk)
%%
%%%%

close all
clear

%% add paths
addpath('tools/')
addpath('tools/VoronoiSphere/')
addpath('tools/TriangleRayIntersection/')
addpath('../API_MO/SOFAtoolbox/')

%% define paths to subject directories (containing recorded sweeps in 'sweeps' folder)
dirlist = [
    "data/20220807-XR-SUBJ001B/";
    ];

%% PROCESSING
for i = 1:length(dirlist)
    subjectdir = convertStringsToChars(dirlist(i));
    
    % get BRIRs from sine sweeps
    inv_pol = false; % change to true if either mics or the loudspeaker has inversed polarisation
    irBank = extractXrIRs(subjectdir,inv_pol);

    % normalize all measurements using max value
    irBank = normalizePeak(irBank);

    % plot all measured impulse responses
    figure('Name','Measured IRs','NumberTitle','off','WindowStyle','docked');
    hold on
    for j = 1:length(irBank)
        plot(irBank(j).fullIR(:,1), '-g')
        if size(irBank(j).fullIR,2) == 2
            plot(irBank(j).fullIR(:,2), '-r')
        end
    end
    
    % save irBank
    save([subjectdir 'irBank.mat'], 'irBank')
    load([subjectdir 'irBank.mat'])

    % create directory for plots
    mkdir([subjectdir 'figures/'])
    plotMagnitudes(irBank, '1-measured', [subjectdir 'figures/'])
    
    % headphone EQ
%     hpEQ(hpirBank, subjectdir)
    
    % HMD influence correction(ITD and magnitude)
%     irBank = hmdCorrection(irBank);

    % time domain windowing
    plotting = 'false';
    irBank = winIRs(irBank, plotting, [subjectdir 'figures/windowing/']); % set 'true' to save plots
    plotMagnitudes(irBank, '2-win', [subjectdir 'figures/'])

    % calculate FF measurement inv filter and normalize all hrirs
    % do the low-frequency extension
    plotting = 'true';
    irBank = normalizeIRs(irBank, plotting, [subjectdir 'figures/']);
    plotMagnitudes(irBank, '3-raw', [subjectdir 'figures/'])

    % diffuse-field equalization
    dfe_enabled = true;
    plotting = 'true';
    irBank = dfeHRIRs(irBank, dfe_enabled, plotting, [subjectdir 'figures/']);
    plotMagnitudes(irBank, '4-dfe', [subjectdir 'figures/'])

    % save sofa file
    saveAsSofa(irBank, subjectdir,'raw')
    saveAsSofa(irBank, subjectdir,'dfe')
    
    % save ambix config file
    addpath('../adt/')
    adt_initialize
    saveAsAmbix(irBank, subjectdir)

    % save processed irBank
    save([subjectdir 'irBankProcessed.mat'], 'irBank')
    load([subjectdir 'irBankProcessed.mat'])
    
    % do barycentric interpolation of hrirs using time alignment
    interpHrirBank = interpHRIRs(irBank,'raw','align');
    saveAsSofa(interpHrirBank,subjectdir,'interp-raw')
    
    interpHrirBank = interpHRIRs(irBank,'dfe','align');
    saveAsSofa(interpHrirBank,subjectdir,'interp-dfe')
end
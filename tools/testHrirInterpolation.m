close all
clear

subjectdir = 'data/20220807-XR-SUBJ001B/';
load([subjectdir 'irBankProcessed.mat'])

% does barycentric interpolation of hrirs
addpath('tools/TriangleRayIntersection/')
addpath('tools/')

interpHrirBank = interpHRIRs(irBank,'raw','minph');
interpHrirBank = interpHRIRs(irBank,'raw','align');

% interpHrirBank = interpHRIRs(irBank,'dfe','minph');
% interpHrirBank = interpHRIRs(irBank,'dfe','align');

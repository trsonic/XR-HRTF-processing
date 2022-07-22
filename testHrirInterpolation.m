close all
clear

subjectdir = 'data/20211126-XR-TR/';
% subjectdir = 'data/20211126-XR-Gavin/';
load([subjectdir 'irBankProcessed.mat'])

% does barycentric interpolation of hrirs
addpath('tools/TriangleRayIntersection/')
addpath('tools/')

interpHrirBank = interpHRIRs(irBank,'raw','minph');
interpHrirBank = interpHRIRs(irBank,'raw','align');

% interpHrirBank = interpHRIRs(irBank,'dfe','minph');
% interpHrirBank = interpHRIRs(irBank,'dfe','align');

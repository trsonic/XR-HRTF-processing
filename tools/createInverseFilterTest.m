close all
clear

subjectdir = 'data/20211126-XR-TR/';
% subjectdir = 'data/20211126-Gavin/';
load([subjectdir 'irBankInvTesting.mat'])

h = irBank(51).winIR(:,1);
Fs = irBank(51).Fs;

% h = [1; zeros(255,1)];
% Fs = 48000;

h = h/max(abs(h));

invh = createInverseFilter(h,Fs);


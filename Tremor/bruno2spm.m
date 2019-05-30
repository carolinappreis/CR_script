ftdata = [];
ftdata.trial{1} = [ONDOP];
ftdata.time{1} = Xaxis;
ftdata.fsample = 1/mean(diff(Xaxis));
ftdata.label = {'M1'};

D = spm_eeg_ft2spm(ftdata,  'ONDOP.mat')
%%
D = chantype(D, 1, 'LFP');

save(D);
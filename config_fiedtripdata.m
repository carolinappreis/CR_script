% data_eeg=struct;


% X = [SNR.bua{1,1}(2:end,:);
fsamp = 500;
X = rand(7,10000);
XT = linspace(0,size(X,2)/fsamp,size(X,2));

data_eeg = [];
data_eeg.label= {'contc1','contc2', 'contc3','contc4','contc5','contc6','contc7'};
data_eeg.trial= {X};
data_eeg.time={XT};
data_eeg.fsample=fsamp;

% 
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [14 20];
filtData = ft_preprocessing(cfg,data_eeg)


cfg = [];
cfg.method = 'hanning';
% cfg.tapsmofrq  = 3;
ft_freqanalysis(cfg,filtData)
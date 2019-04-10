% dat{1} = randn(10, 100);
% dat{2} = randn(10, 100);
% 
% dat{2}(:, 40:60) = dat{2}(:, 40:60)+2;
% 
% fsample = 250;

dat{1}=A;
dat{2}=B;%zeros(size(A,1),size(A,2));
fsample=1000;

for i = 1:numel(dat)
    timelock(i).dimord = 'rpt_chan_time';
    timelock(i).label = {'Ch1'};
    timelock(i).time = (1:size(dat{i}, 2))/fsample;
    timelock(i).trial = reshape(dat{i}, [size(dat{i}, 1), 1, size(dat{i}, 2)]);
end

cfg=[];
cfg.latency   = 'all';
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.tail = 0;
cfg.alpha = 0.05;
cfg.numrandomization = 2000;

cfg.correctm = 'cluster'; %'no|max|cluster|bonferoni|holms|fdr'

cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';

design = zeros(1, size(dat{1}, 1)+size(dat{2}, 1));
design(1,1:size(dat{1}, 1)) = 1;
design(1,(size(dat{1}, 1)+1):end)=2;
design(2, :) = [1:size(dat{1}, 1) 1:size(dat{2}, 1)];


cfg.design   = design;   % design matrix
cfg.ivar  = 1;  % number or list with indices, independent variable(s)
cfg.uvar  = 2;

cfg.neighbours=[];
for c=1:length(timelock(1).label)
    cfg.neighbours(c).label = timelock(1).label{c};
    cfg.neighbours(c).neighblabel ={};
end

cfg.parameter='trial';
%%
stats= ft_timelockstatistics(cfg, timelock(1), timelock(2));

% stats.posclusters(1)
% stats.negclusters(1)
% 
% stats.posclusterslabelmat
% stats.negclusterslabelmat
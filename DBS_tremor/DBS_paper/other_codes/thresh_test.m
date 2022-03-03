clearvars -except envelope_all
clear all
close all
cond={'NS';'HF';'C'};
cohort=[ 1 3 4 6];

iii=2;
trial=2;
co=1;

load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/clean_SH_spirals/P0',num2str(cohort(iii)),'_clean_',num2str(cond{co,1}),num2str(trial),'_SH.mat'));


t_spi=cell(size(cohort,2),2,size(cond,1));

% t_spi{2,1,1}=[[1 1939];[1960 4304];];
t_spi{2,1,1}=[1 length(tempo)];
t_spi{2,2,1}=[[1 2534];[2534 5228];[5228 length(tempo)]];
t_spi{2,1,2}=[[1 1892];[1894 length(tempo)]];
t_spi{2,1,3}=[[47 2105];[2105 4072];[4072 6236]];


t_spi{3,1,1}=[[1 1687];[1687 2947];[2947 4354];[4354 5601];[5650 length(tempo)]];
t_spi{3,2,1}=[[1 1668];[1668 2979];[2979 4341];[4341 length(tempo)]];
t_spi{3,1,2}=[[1 1130];[1130 2517];[2517 4207];[4207 5605];[5605 length(tempo)]];
t_spi{3,1,3}=[[1 1420];[1420 2647];[2647 3834];[3834 4923];[4923 6132];[6132 7283];[7283 8507];[8507 9827];[9827 10990];[10990 length(tempo)]];

t_spi{4,1,1}=[[1 6260];[7209 length(tempo)]];
t_spi{4,1,2}=[[1 5236];[5236 length(tempo)]];
t_spi{4,1,3}=[[1 5932];[5932 length(tempo)]];



%%% filt & var
[a,b]=  butter(2, [2/(0.5*samplerate2) 8/(0.5*samplerate2)], 'bandpass'); %15
filt_cs=(filtfilt(a,b,signal(3,:)));
env=abs(hilbert(filt_cs));

senv=smoothdata(env,'movmean',samplerate2);

% plot3(tempo,signal(1,:),signal(2,:),'.')

%%% spatial
% signal=data;
centre=[535 361];
d1=signal(1,:)-centre(1);
d2=signal(2,:)-centre(2);
[px,py]=cart2pol(d1,d2);
polarscatter(px(:), py(:), [], senv(:),'filled')
colorbar
colormap winter
title('envelope')


figure(5)
for mm=1:size(t_spi{iii,trial,co}(:,1),1)
    
% subplot(1,size(t_spi{iii,trial,co}(:,1),1),mm)
%      subplot(2,5,mm)
  subplot(2,3,mm)
%    subplot(1,3,mm)
    polarscatter(px(t_spi{iii,trial,co}(mm,1):t_spi{iii,trial,co}(mm,2)),py(t_spi{iii,trial,co}(mm,1):t_spi{iii,trial,co}(mm,2)),[],senv(t_spi{iii,trial,co}(mm,1):t_spi{iii,trial,co}(mm,2)),'filled');
    colorbar
    colormap winter
%     load (strcat('env_iii0',num2str((iii)))); caxis([min(envelope_all)/mean(envelope_all) max(envelope_all)/mean(envelope_all)])
       caxis([min(env) max(env)-5])
%       caxis([min(env(t_spi{iii,trial,co}(mm,1):t_spi{iii,trial,co}(mm,2))) max(env(t_spi{iii,trial,co}(mm,1):t_spi{iii,trial,co}(mm,2)))])
%   caxis([min(env(t_spi{iii,trial,co}(mm,1):t_spi{iii,trial,co}(mm,2)))./median(env) max(env(t_spi{iii,trial,co}(mm,1):t_spi{iii,trial,co}(mm,2)))./median(env)])
end

%     figure
%     plot(envelope_all)
%     box('off')
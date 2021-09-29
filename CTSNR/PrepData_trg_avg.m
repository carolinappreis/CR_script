%% JUXTA
clear all; close all
 cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/Juxta SUA_act_mat/mat')
load ('SUA_BZ')

data=data_all;
ecog=WaveData_DC; 

%% SUA
clear all; close all
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/probe SUA_act_mat')
load('data_SUA_SNR.mat')

data=cell2mat(data_region);
ecog=[];
for i=1:size(data_region,1)
    ecog = [ecog ; repmat(Ecog_region(i,:),size(data_region{i,1},1),1)];
end

%% Analyses
clearvars -except data ecog name
cd('/Users/Carolina/Documents/GitHub/CR_script/CTSNR')

[fig]=trg_avg(data,ecog,name);
[pref_in, pref_out, mean_ang,spike_rate]=pref_lock_ang(data,ecog,name)
[fig]=broad_cycy(data,ecog,name);
[fig]=lock_prop(data,ecog,name);

%%%% extras

% check firing rate units
% % [srate]=firing_rate(data_region);

%%% envelope variability; burst metrics
% srn=1000;
% [b,a]=butter(2,[15/(0.5*srn) 35/(0.5*srn)],'bandpass');
% for i=1:size(Ecog_region,1)
%     ctx(i,:)=filtfilt(b,a,Ecog_region(i,:));
% end
% [env_var]=burst_var(ctx);
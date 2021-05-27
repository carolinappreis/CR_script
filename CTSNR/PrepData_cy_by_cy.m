%% JUXTA
clear all; close all
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/Juxta SUA_act_mat/mat')
load('BZ_cycle')

data=units_match;
ecog=ecogbf_match;

%% SUA
clear all; close all
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\probe SUA_act_mat')
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/probe SUA_act_mat')
load ('NEW_SNR_cycle.mat')

data_region=units_match(~cellfun('isempty',units_match));
Ecog_region=ecogbf_match(any(ecogbf_match,2),:);

data=cell2mat(data_region);
ecog=[];
for i=1:size(data_region,1)
    ecog = [ecog ; repmat(Ecog_region(i,:),size(data_region{i,1},1),1)];
end

%% Analyses

clearvars -except data ecog name
[fig]=cy_by_cy(data,ecog,name);
[fig]=lock_prop(data,ecog,name);


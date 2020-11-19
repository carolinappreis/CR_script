clear all
close all
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/final_mats')
load('BZ_SNR.mat');

for i=1:size(rat_labels,1)
    Ecog{i,1}=BZ_bua{i,1}(1,:);
end


refs={'Ecog' 'BZ_bua' 'SN_bua'};




clear all
close all
% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\cleaned_rc12_noaddon.mat')
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cleaned_rc12_noaddon.mat')

for i=1:size(tt1,1)
    for ax =1:size(tt1,2)

   abs_change(i,ax,:)=(abs(nanmedian(tt1{i,ax})));
   sum_change(i,ax)=sum(abs(nanmedian(tt1{i,ax})));
    end
end
   


[p,tbl,stats] = anova1((squeeze(abs_change(3,:,:))'))
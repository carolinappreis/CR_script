
clear all

load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_final_spiral.mat');
A=[skewness(out.mod_amp{4,2}{2,1}(:)) skewness(out.mod_amp{3,2}{3,1}(:)) skewness(out.mod_amp{2,2}{3,1}(:)) skewness(out.mod_amp{1,2}{3,1}(:))];
B=[skewness(out.mod_amp{4,1}(2,:)) skewness(out.mod_amp{3,1}(3,:)) skewness(out.mod_amp{2,1}(3,:)) skewness(out.mod_amp{1,1}(3,:))];


A=[(out.mod_amp{4,2}{2,1}(:)); (out.mod_amp{3,2}{3,1}(:)); (out.mod_amp{2,2}{3,1}(:)) ;(out.mod_amp{1,2}{3,1}(:))];
B=[(out.mod_amp{4,1}(2,:)) (out.mod_amp{3,1}(3,:)) (out.mod_amp{2,1}(3,:)) (out.mod_amp{1,1}(3,:))];


[h,p]=ttest(A,B);


clear all
load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_final_posture.mat');
A=[skewness(out.mod_amp{4,2}{1,1}(:)) skewness(out.mod_amp{3,2}{3,1}(:)) skewness(out.mod_amp{2,2}{3,1}(:)) skewness(out.mod_amp{1,2}{3,1}(:))];
B=[skewness(out.mod_amp{4,1}(1,:)) skewness(out.mod_amp{3,1}(3,:)) skewness(out.mod_amp{2,1}(3,:)) skewness(out.mod_amp{1,1}(3,:))];


clear all
load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_final_posture.mat');
A=[(out.mod_amp{4,2}{1,1}(:)) (out.mod_amp{3,2}{3,1}(:)) (out.mod_amp{2,2}{3,1}(:)) (out.mod_amp{1,2}{3,1}(:))];
B=[(out.mod_amp{4,1}(1,:)) (out.mod_amp{3,1}(3,:)) (out.mod_amp{2,1}(3,:)) (out.mod_amp{1,1}(3,:))];


[h,p]=ttest(A,B)




clear all
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cluster_out_mc.mat');
A=[]; B=[];
for i=1:size(out.change_c,1)
    A=[A skewness(out.change_c{i,2}{1,1}(:))'];
    B=[B skewness(out.change_c{i,1}(1,:))'];
end

[h,p]=ttest(B,A)


clear all
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cluster_out_mc.mat');
A=[]; B=[];
for i=1:size(out.change_c,1)
    A=[A (out.change_c{i,2}{1,1}(:))'];
    B=[B (out.change_c{i,1}(1,:))'];
end


subplot(1,2,1)
histogram(A)
xlim([-1 2])
xlabel('change in tremor severity')
ylabel('counts')
box('off')
subplot(1,2,2)
histogram(B)
xlim([-1 2])
xlabel('change in tremor severity')
ylabel('counts')
box('off')
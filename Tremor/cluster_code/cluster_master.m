clear all; clear all
close all
cohort = [ 2 3 4 5 8 10 11 13 16 17];

% nostim = NaN(length(cohort),3,1e6);
cond={'NS';'RS';'PLS'};
clust=struct;
start=cell(10,3);
ending=cell(10,3);
xx=cell(10,1);

for iii = 1
%     1:length(cohort)
    clearvars -except  cohort cond iii clust s start ending
    all=[];
    for co=1:size(cond,1)
        
        load(strcat('/Users/Carolina/OneDrive - Nexus365/PERI-STIM/DATA/p0',num2str(cohort(iii)),'_',num2str(cond{co,1}),'.mat'))
        
        [d]=preprocess(SmrData);
        
        samplerateold=d.samplerateold; samplerate=d.samplerate;

        data=d.data_raw;  tremor_ds=d.data_ds;
        [afilt, bfilt, start, ending, xx ]= startend(data,tremor_ds,samplerateold,samplerate,iii,co,start,ending); clear data tremor_ds samplerateold;
        
        data=d.data_ds;
        [s]=zfiltenv(data,bfilt,afilt,co,iii); clear data bfilt afilt;
        

        filt_tremor=s.filt;
        [pc_trials,runs,start,ending]=pcomp(filt_tremor,start,ending,co,samplerate,iii);
        
        all= [all ; pc_trials]; clear pc seg_pc d filt_tremor
        
        pc_comp{iii,co}=pc_trials;
    end
    
    [clust]=clustering(all,iii,runs,clust,pc_comp);
    
    [tt,segments]=modulation(clust,s);
    
end
clearvars -except clust
save('cluster_master_out.mat')

% 
% load('/Users/Carolina/Documents/GitHub/CR_script/Tremor/cluster_code/cluster_master_out.mat')
% 
% for iii=1:10

% end
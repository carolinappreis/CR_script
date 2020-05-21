clear all; close all
cohort = [ 2 3 4 5 8 10 11 13 16 17];
cond={'NS';'RS';'PLS'};
clust=struct; out=struct; start=cell(10,3); ending=cell(10,3); yy=cell(10,3);

for iii = 1
    %     :length(cohort)
    clearvars -except  cohort cond iii clust s start ending yy out
    all=[];
    for co=1:size(cond,1)
        
        load(strcat('/Users/Carolina/OneDrive - Nexus365/PERI-STIM/DATA/p0',num2str(cohort(iii)),'_',num2str(cond{co,1}),'.mat'))
        
        [d]=preprocess(SmrData);
        
        samplerateold=d.samplerateold; samplerate=d.samplerate;
        
        data=d.data_raw;  tremor_ds=d.data_ds;
        [afilt, bfilt, start, ending, yy]= startend(data,tremor_ds,samplerateold,samplerate,iii,co,start,ending,yy); clear data tremor_ds samplerateold;
        
        data=d.data_ds;
        [s]=zfiltenv(data,bfilt,afilt,co,iii); clear data bfilt afilt;
        
        [pc_trials,start,ending,out,yy]=mod_nc(start,ending,co,samplerate,iii,s,yy,out);
        
        [nc]=plots_nc(out,iii);
        
        runs(1,co)=length(start{iii,co});
        all= [all ; pc_trials];
        pc_comp{iii,co}=pc_trials;
    end
    
%     [clust]=clustering(all,iii,runs,clust,pc_comp);
    [clust]=clustering2(all,iii,runs,clust,pc_comp);
    
    for co=1:3
        [out,s]=mod_c(clust,s,samplerate,out,yy,start,ending,co,iii);
    end
    
    [c1]=plots_c(out,iii);
end

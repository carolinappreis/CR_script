clear; close 
cohort = [ 2 3 4 5 8 10 11 13 16 17];
cond={'NS';'RS';'PLS'};
clust=struct; out=struct; start=cell(10,3); ending=cell(10,3); yy=cell(10,3); h_up=cell(10,3); s=struct;
rng('default')
gen=(rng);
% 
%   load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cluster_trials.mat','clust_trials');clust.idx=clust_trials;
%  load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cluster_master_output.mat');
for iii =1:length(cohort)
    clearvars -except  cohort cond iii clust s start ending yy out gen h_up
    all=[];
    for co=1:size(cond,1)
        
        load(strcat('/Users/Carolina/OneDrive - Nexus365/PERI-STIM/DATA/P0',num2str(cohort(iii)),'_',num2str(cond{co,1}),'.mat'))
        
        [d]=preprocess(SmrData); samplerateold=d.samplerateold; samplerate=d.samplerate;
        
%         [h]=twitch(iii,d,samplerateold,samplerate,co); clear h;
        
        [afilt, bfilt, start, ending, yy, h_up]= startend(d,samplerateold,samplerate,iii,co,start,ending,yy,h_up);
        
        [s]=zfiltenv(d,bfilt,afilt,co,iii,s); clear afilt bfilt
        
        % [start,ending,out,yy]=mod_nc(start,ending,co,samplerate,iii,s,yy,out,gen); %%%tremor amplitude change without clustering
        
    end
    
   [clust,out]=clustering2(out,iii,clust,start,ending,yy); %% clustering analysis
    
    for co=1:size(cond,1)
        [out]=mod_c(clust,out,co,iii,s,h_up);   %%%% tremor amplitude change with clustering
    end 
end
clearvars -except out

[out]=polyfit(out); %non uniformity of arc's polynomial fits with r2 and f-stats
[r,p]=plots_c(out);  %plots ARCs with clustering
[out]=psd_rs(out);  %power of RS segments 
stats=struct;
[stats]=pwelch3(out,stats);  % stat comparison between power and freq between amplification supression and no stim
[stats]=group_aligned(out,stats);  % group ARC - stats 
% % [nc]=plots_nc(out); %plots ARCs without clustering ---- needs fx mod_nc






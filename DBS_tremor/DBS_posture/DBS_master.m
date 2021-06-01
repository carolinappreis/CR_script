clear; close
cohort = [ 1 3 4 6];
% cond={'NS';'RS';'PLS'};
cond={'NS';'RS'};
clust=struct; out=struct; start=cell(10,3); ending=cell(10,3); yy=cell(10,3); h_up=cell(10,3); s=struct;
rng('default')
gen=(rng);

spiral=0;
 load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_cluster_out.mat');

for iii =  1:length(cohort)
    clearvars -except  cohort cond iii clust s start ending yy out gen h_up spiral
    for co=1:size(cond,1)
        
        load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA/P0',num2str(cohort(iii)),'_',num2str(cond{co,1}),'.mat'))
        
        [d]=dbs_preprocess(SmrData); samplerateold=d.samplerateold; samplerate=d.samplerate;
        
        [peak_ax, start, ending, yy, h_up]= dbs_startend(d,samplerateold,samplerate,iii,co,start,ending,yy,h_up,spiral);
        
        [s]=dbs_zfiltenv(d,peak_ax,co,iii,s,samplerate);
        
%                 if ~isnan(clust.win(iii,1))
%                     [out]=cluster_intime(clust,s,iii,co,out);
%                 end
%         
%          [start,ending,out,yy]=dbs_mod_nc(start,ending,co,samplerate,iii,s,yy,out); %%%tremor amplitude change without clustering
    end
    
  [clust,out]=dbs_clustering2(out,iii,clust,start,ending,yy,spiral); %% clustering analysis
%     
 for co=1:size(cond,1)
     [out]=mod_c(clust,out,co,iii,s,h_up,cohort,spiral);   %%%% tremor amplitude change with clustering
 end
    
end
%   [nc]=dbs_plots_nc(out)
  clearvars -except out clust

% clear
% load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_cluster_out.mat','out');


[out]=polyfit(out); %non uniformity of arc's polynomial fits with r2 and f-stats
[r,p]=plots_c(out);  %plots ARCs with clustering

stats=struct;
[out]=psd_rs(out);  %power of RS segments 
[stats]=pwelch3(out,stats);  % have to run psd_rs first % stat comparison between power and freq between amplification supression and no stim

[stats]=group_aligned(out,stats);  % group ARC - stats
% % [nc]=plots_nc(out); %plots ARCs without clustering ---- needs fx mod_nc


%%% continuos stim during posture
% continuous_posture (change with hand up, tap, hand down...)
% continuos_p_pca (pca contribution and main axis at each stage of the task
% power_pls_hfs_posture (power spectra hfs vs continuous stim, 2pts)
%ns_rs_pls_hf_posture (psd all conditions for 2 pts)
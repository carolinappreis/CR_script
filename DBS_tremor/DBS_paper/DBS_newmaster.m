cd('/Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/DBS_paper')
clear; close
cohort = [ 1 3 4 6];
cond={'NS';'RS'};
clust=struct; out=struct; start=cell(10,1); ending=cell(10,1); yy=cell(10,3); h_up=cell(10,3); s=struct;
rng('default')
gen=(rng);

spiral=0;
% if spiral==0
%    load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_cluster_out.mat');
% else
%     load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_clust_out_sp.mat');
% end

for iii = 1:length(cohort)
    clearvars -except  cohort cond iii clust s start ending yy out gen h_up spiral
    for co=1:size(cond,1)
        
        load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA/P0',num2str(cohort(iii)),'_',num2str(cond{co,1}),'.mat'))
        
        [d]=dbs_preprocess(SmrData); samplerateold=d.samplerateold; samplerate=d.samplerate;
        
        [peak_ax, start, ending, yy, h_up]= srtend(d,samplerateold,samplerate,iii,co,start,ending,yy,h_up,spiral);
        
        [s]=zfiltenv_simple(d,peak_ax,co,iii,s,samplerate);
        
        [out]=mod_new(clust,out,co,iii,s,spiral,start,ending);  %%%% tremor amplitude change with clustering
    end
    
    [clust,out]=clust_inside(out,iii,clust,start,ending,yy,spiral);
    
    
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
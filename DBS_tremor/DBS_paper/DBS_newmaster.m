cd('/Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/DBS_paper')
clear; close
cohort = [ 1 3 4 6];
cond={'NS';'RS'};
clust=struct; out=struct; start=cell(10,1); ending=cell(10,1); yy=cell(10,3); h_up=cell(10,3); s=struct;
rng('default')
gen=(rng);

spiral=0;

% if spiral==0
%   load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_cluster_out.mat');
% else
%     load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_clust_out_sp.mat');
% end

for iii = 1
%     :length(cohort)
    clearvars -except  cohort cond iii clust s start ending yy out gen h_up spiral
    
    for co=1
%         :size(cond,1)  
    load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA/P0',num2str(cohort(iii)),'_',num2str(cond{co,1}),'.mat'))
    [d]=dbs_preprocess(SmrData); samplerateold=d.samplerateold; samplerate=d.samplerate;
    [peak_ax, start, ending, yy, h_up]= srtend(d,samplerateold,samplerate,iii,co,start,ending,yy,h_up,spiral);
    [s]=zfiltenv_simple(d,peak_ax,co,iii,s,samplerate);
%   [freq_bl,amp_bl]=pre_mod(co,iii,s,spiral); 
    end
    
    load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/a1')
    load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/f1')


    [clust,out]=clust_inside(out,iii,clust,spiral);
    
    for co=1:size(cond,1)
    load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA/P0',num2str(cohort(iii)),'_',num2str(cond{co,1}),'.mat'))
    [d]=dbs_preprocess(SmrData); samplerateold=d.samplerateold; samplerate=d.samplerate;
    [peak_ax, start, ending, yy, h_up]= srtend(d,samplerateold,samplerate,iii,co,start,ending,yy,h_up,spiral);
    [s]=zfiltenv_simple(d,peak_ax,co,iii,s,samplerate);
    [out,clust]=mod2(out,co,iii,s,freq_bl,amp_bl,start,ending,yy,clust,spiral);
    end
end
clearvars -except out clust


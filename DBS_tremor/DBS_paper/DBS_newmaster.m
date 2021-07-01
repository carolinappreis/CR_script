cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA')
cd('/Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/DBS_paper')
clear; close
cohort = [ 1 3 4 6];
cond={'NS';'RS'};
clust=struct; out=struct; start=cell(10,1); ending=cell(10,1); yy=cell(10,3); h_up=cell(10,3); s=struct; freq_bl=[]; amp_bl=[];
rng('default')
gen=(rng);
spiral=0;

if spiral==0
    load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_final_posture.mat');
else
    load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_final_spiral.mat');
end


for iii = 2:length(cohort)
    clearvars -except  cohort cond iii clust s start ending yy out gen h_up spiral freq_bl amp_bl
    
    % %     for co=1:size(cond,1)
    % %     load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA/P0',num2str(cohort(iii)),'_',num2str(cond{co,1}),'.mat'))
    % %     [d]=dbs_preprocess(SmrData); samplerateold=d.samplerateold; samplerate=d.samplerate;
    % %     [peak_ax, start, ending, yy, h_up]= srtend(d,samplerateold,samplerate,iii,co,start,ending,yy,h_up,spiral);
    % %     [s]=zfiltenv_simple(d,peak_ax,co,iii,s,samplerate);
    % %     [freq_bl,amp_bl,out]=pre_mod(co,iii,s,spiral,start,ending,out,freq_bl,amp_bl);
    % %     end
    % %
    % %     [clust,out]=clust_inside(out,iii,clust,spiral);
    
    for co=1
%         :size(cond,1)
        load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA/P0',num2str(cohort(iii)),'_',num2str(cond{co,1}),'.mat'))
        [d]=dbs_preprocess(SmrData); samplerateold=d.samplerateold; samplerate=d.samplerate;
        [peak_ax, start, ending, yy, h_up]= srtend(d,samplerateold,samplerate,iii,co,start,ending,yy,h_up,spiral);
        [s]=zfiltenv_simple(d,peak_ax,co,iii,s,samplerate);
        [out,clust,freq_bl,amp_bl]=mod2(out,co,iii,s,freq_bl,amp_bl,start,ending,yy,clust,spiral,h_up);
    end
end
clearvars -except out clust

% if spiral==0
%   load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_final_posture.mat');
% else
%   load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_final_spiral.mat');
% end


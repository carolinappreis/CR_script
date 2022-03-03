cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA')
cd('/Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/DBS_paper')
clear; close
cohort = [ 1 3 4 6];
cond={'NS';'RS'};
clust=struct; out=struct; start=cell(10,1); ending=cell(10,1); yy=cell(10,3); h_up=cell(10,3); s=struct; freq_bl=[]; amp_bl=[];
rng('default')
gen=(rng);
spiral=1;

if spiral==0
    load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_final_posture_unclust.mat');
else
    load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_final_spiral_unclust.mat');
end


for iii =  1:length(cohort)
    clearvars -except  cohort cond iii clust s start ending yy out gen h_up spiral freq_bl amp_bl avg power dat
    
        for co=1:size(cond,1)
        load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA/P0',num2str(cohort(iii)),'_',num2str(cond{co,1}),'.mat'))
        [d]=dbs_preprocess(SmrData); samplerateold=d.samplerateold; samplerate=d.samplerate;
        [peak_ax, start, ending, yy, h_up]= srtend(d,samplerateold,samplerate,iii,co,start,ending,yy,h_up,spiral);
        [s]=zfiltenv_simple(d,peak_ax,co,iii,s,samplerate);
          
         [match_ax]=link_ax(spiral);
%          [Pxx,F] = pwelch(s.raw{iii, 1} (match_ax(co,iii,1),h_up{iii,1}), samplerate, [], samplerate, samplerate);
%           dat{iii,1}=s.raw{iii, 1} (match_ax(co,iii,1),h_up{iii,1});       
%             power(iii,:)=Pxx; clear Pxx
         
         
         avg(spiral+1,iii,:)=[mean(s.env_acc{iii,co}(match_ax(1,iii,1),h_up{iii,co}),2) std(s.env_acc{iii,co}(match_ax(1,iii,1),h_up{iii,co}))];
         
         [freq_bl,amp_bl,out]=pre_mod(co,iii,s,spiral,start,ending,out,freq_bl,amp_bl);
        end
    
%          [clust,out]=clust_inside(out,iii,clust,spiral);
    
      for co=1:size(cond,1)
        load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA/P0',num2str(cohort(iii)),'_',num2str(cond{co,1}),'.mat'))
        [d]=dbs_preprocess(SmrData); samplerateold=d.samplerateold; samplerate=d.samplerate;
        [peak_ax, start, ending, yy, h_up]= srtend(d,samplerateold,samplerate,iii,co,start,ending,yy,h_up,spiral);
        [s]=zfiltenv_simple(d,peak_ax,co,iii,s,samplerate);
        [out,clust,freq_bl,amp_bl]=mod2(out,co,iii,s,freq_bl,amp_bl,start,ending,yy,clust,spiral,h_up);
     end
end
clearvars -except out clust spiral

cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA')

% if spiral==0
%   load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_final_posture.mat');
% else
%   load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_final_spiral.mat');
% end

%%%% code loop through the 10 patients and 2 conditions (no-stim and random
%%%% stim (blocks of stim during 5 seconds at 12 different phases)

clear; close
cohort = [ 2 3 4 5 8 10 11 13 16 17];
cond={'NS';'RS';'PLS'};
clust=struct; out=struct; start=cell(10,3); ending=cell(10,3); yy=cell(10,3); h_up=cell(10,3); s=struct;
% rng('default')
% gen=(rng);
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cluster_out_mc.mat');

for iii =1: length(cohort)
    clearvars -except  cohort cond iii clust s start ending yy out gen h_up
    
   for co= 1 % NS and RS conditions       
        load(strcat('/Users/Carolina/OneDrive - Nexus365/PERI-STIM/DATA/P0',num2str(cohort(iii)),'_',num2str(cond{co,1}),'.mat'))

        [d]=preprocess(SmrData); samplerateold=d.samplerateold; samplerate=d.samplerate;
        
%          [h]=twitch(iii,d,samplerateold,samplerate,co); clear h; %% jolt induced by stim recorded on main axis of the 3axial acc and emg for patients 9 and 10
        
        [afilt, bfilt, start, ending, yy, h_up]= startend(d,samplerateold,samplerate,iii,co,start,ending,yy,h_up); %% start and end points of stim runs
        
        [s]=zfiltenv(d,bfilt,afilt,co,iii,s); clear afilt bfilt %% zscore / filter / envelope / phase
    end
    
%         [clust,out]=clustering2(out,iii,clust,start,ending,yy); %% clustering analyses
    
     for co=1
        [out]=mod_c(clust,out,co,iii,s,cohort,h_up);   %%%% tremor modulation
    end
   
    
end
clearvars -except out clust
cd('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data')
%%%%% save('cluster_out_mc.mat');

%
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cluster_out_mc.mat');
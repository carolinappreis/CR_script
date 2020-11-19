
clear all
close all
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/final_mats')
load 'animal_lesion_nolesion.mat' %%% from code sort_lesion_nonlesion.m (microsoft version needs adapting for mac)
state_mat=lesion;% select subject group lesion vs no lesioned
nfile_i=3:(newfile);
med_thal_cell=cell(1,1);
med_ctx=cell(1,1);
n_bursts=cell(1,1);
x=0;
check=[];
for nn=1:length(state_mat)
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/final_mats/KN')
    clearvars -except x med_thal_cell med_ctx nfile_i state_mat nn A lesion data rat_label name
    label=A(state_mat(nn),1:(find(A((state_mat((nn))),:)=='.')-1));
    load(label)
    B=who;

    %identify EEG signal and Probe signal
    ctxchan=[];
    thalchan=[];
    for ii=1:size(B,1)
        if ~isempty(min(find(B{ii}=='i')) & min(find(B{ii}=='p')) & min(find(B{ii}=='s')))
            ctxchan=ii;
        elseif ~isempty(min(find(B{ii}=='p')) & min(find(B{ii}=='r')) & min(find(B{ii}=='o')) & min(find(B{ii}=='b')))
            thalchan=[thalchan ii];
        end
    end
    
    %Discriminate probe region
    thal_BZ=[];%basal-ganglia recipient thalamus
    thal_CZ=[];%cerebellar recipient thalamus
    thal_ZI=[];%Zona incerta
    thal_AV=[];%Amteroventral thalamus
    thal_SNr=[];%Substantia Nigra reticulata
    thal_NW=[];%No where
    thal_PC=[];%Para-central thalamus
    thal_AM=[];%Antero medial thalamus
    thal_Rt=[];%Retucular thalamic nucleus
    for i=1:length(thalchan)
        eval(['location=' B{thalchan(i)} '.location']);
        if ~isempty(min(find(location=='B')) & min(find(location=='Z')));
            thal_BZ=[thal_BZ,i];
        elseif ~isempty(min(find(location=='C')) & min(find(location=='Z')));
            thal_CZ=[thal_CZ,i];
        elseif ~isempty(min(find(location=='A')) & min(find(location=='V')));
            thal_AV=[thal_AV,i];
        elseif ~isempty(min(find(location=='Z')) & min(find(location=='I')));
            thal_ZI=[thal_ZI,i];
        elseif ~isempty(min(find(location=='S')) & min(find(location=='N')) & min(find(location=='r')));
            thal_SNr=[thal_SNr,i];
        elseif ~isempty(min(find(location=='n')) & min(find(location=='o')) & min(find(location=='w')));
            thal_NW=[thal_NW,i];
        elseif ~isempty(min(find(location=='P')) & min(find(location=='C')));
            thal_PC=[thal_PC,i];
        elseif ~isempty(min(find(location=='A')) & min(find(location=='M')));
            thal_AM=[thal_AM,i];
        elseif ~isempty(min(find(location=='R')) & min(find(location=='t')));
            thal_Rt=[thal_Rt,i];
        end
    end
    
    choose=thal_SNr;
        
    if ~isempty(choose)    
    x=x+1;
    chanofinterest= thalchan(choose); %%% Chose thal location
    

    eval(['samprateold=1/' B{ctxchan} '.interval;']);
    eval(['WaveData(1,:)=' B{ctxchan} '.values;']);
    WaveData=double(WaveData);
    
    
    ts=timeseries(WaveData(1,:),0:(1/samprateold):((size(WaveData,2)-1)/samprateold));
    ts1=resample(ts,0:0.001:((size(WaveData,2)-1)/samprateold),'linear');
    WaveData_DC(1,:)=ts1.data;
    timeold=0:(1/samprateold):(length(WaveData(1,:))-1)*(1/samprateold);
    for ii=1:length(chanofinterest)
        eval(['WaveData(1+ii,:)=' B{chanofinterest(ii)} '.values;']);
        WaveData=double(WaveData);
        WaveData_DC(1+ii,:)=makemua_CR_1(WaveData(1+ii,:)',0.001,0.003,round(samprateold),1000,4);
    end
   
    bua{x,1}=WaveData_DC;
    rat_label{x,1}=A(state_mat(nn),:);
    end
end







%a)select the region we want to analyse
%b)do the spike sorting - BUA
%c)filter the EGG in the beta range
%d)leave thal region unfiltered

clear all
cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A2_Thal/MAT')
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A2_Thal\mat')

load 'animal_lesion_nolesion.mat'

state_mat=lesion;% lesion vs nolesioned
nfile_i=3:(newfile);
rr=1;
for nn=1:length(state_mat);
    
    cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/MATLAB/KN')
    %  cd ('\Users\creis\Documents\MATLAB\KN')
    clearvars -except nfile_i state_mat A rr MUA ctx nn ctx_MUA animals lesion
    name=A(state_mat(nn),1:(find(A((state_mat((nn))),:)=='.')-1))
    load(name)
    B=who;
    
    ctxchan=[];
    thalchan=[];
    for ii=1:size(B,1)
        if ~isempty(min(find(B{ii}=='i')) & min(find(B{ii}=='p')) & min(find(B{ii}=='s')))
            ctxchan=ii;
        elseif ~isempty(min(find(B{ii}=='p')) & min(find(B{ii}=='r')) & min(find(B{ii}=='o')) & min(find(B{ii}=='b')))
            thalchan=[thalchan ii];
        end
    end
    
    thal_BZ=[];%basal-ganglia recipient thalamus
    thal_SNr=[];%Substantia Nigra reticulata
    thal_CZ=[];
    
    for i=1:length(thalchan)
        eval(['location=' B{thalchan(i)} '.location']);
        if ~isempty(min(find(location=='B')) & min(find(location=='Z')));
            thal_BZ=[thal_BZ,i];
        elseif ~isempty(min(find(location=='S')) & min(find(location=='N')) & min(find(location=='r')));
            thal_SNr=[thal_SNr,i];
        elseif ~isempty(min(find(location=='C')) & min(find(location=='Z')));
            thal_CZ=[thal_CZ,i];
        end
    end
    
    chanofinterest= thalchan(thal_SNr);
    eval(['samprateold=1/' B{ctxchan} '.interval;']);
    eval(['WaveData(1,:)=' B{ctxchan} '.values;']);
    WaveData=double(WaveData);
    
    cd ('/Users/Carolina/Documents/GitHub/CR_script/SUA/probe SUA_act_mat')
    
    %  cd ('C:\Users\creis\Documents\GitHub\CR_script\SUA\probe SUA_act_mat')
    
    ts=timeseries(WaveData(1,:),0:(1/samprateold):((size(WaveData,2)-1)/samprateold));
    ts1=resample(ts,0:0.001:((size(WaveData,2)-1)/samprateold),'linear');
    ctx_MUA(1,:)=ts1.data;
    timeold=0:(1/samprateold):(length(WaveData(1,:))-1)*(1/samprateold);
    for ii=1:length(chanofinterest)
        eval(['WaveData(1+ii,:)=' B{chanofinterest(ii)} '.values;']);
        WaveData=double(WaveData);
        ctx_MUA(1+ii,:)=makemua_CR_spikes(WaveData(1+ii,:)',0.00005,0.00005,round(samprateold),1000,5);
    end
    MUA{nn,1}=ctx_MUA(2:(end),:);
    ctx(nn,:)=ctx_MUA(1,:);
    clear ctx_MUA
end

n=[];
for i=1:size(MUA,1)
    if ~isempty (MUA{i,1})
        n=[n i];
    end
end

animals=A(lesion(n),:);
% 
% clearvars -except MUA ctx animals
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\MUA')
% save 'newMUA_SNR.mat'




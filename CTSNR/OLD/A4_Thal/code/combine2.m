clear all
  cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
% cd ('\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
load ('animal_lesion_nolesion.mat','A','lesion')
load 'CZ_ctx_probe.mat'


n=[];
for i=1:29
if size(WaveData_DCall{i,1},1)>1
n=[n i];
end
end


clear all
close all
  cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
load 'CZ_ctx_probe'
fq=15:35;

samprate=1000;
ctx_b1=[];
clust_b=[];
clust_s=[];
ctcts=0;
for r=1:size(data,1); %electrodes
%     thal_contact=[];
    coh_thal=[];
    filt_thal1=[];
    power_coh1=[];
    phase_thal1=[];
    env_thal1=[];
    power=[];
    coherence1=[];
    bua_coh1=[];
    m=1;
    [power1,f1]=pwelch(data{r,1}(1,:),1000,[],1000,1000);
    for rr=2:size(data{r,1},1) %contacts
        [power,f]=pwelch(data{r,1}(rr,:),1000,[],1000,1000); %ctx power spectracontact power spectra
        [Pxx_ind,F_ind]=mscohere(data{r,1}(1,:),data{r,1}(rr,:),samprate,[],samprate,samprate); %Magnitude-squared coherence between ctx-a given contact
       sum(Pxx_ind(fq))./sum(Pxx_ind(1:end))
        if  sum(Pxx_ind(fq))./sum(Pxx_ind(1:end))>0.1
            freq= find(Pxx_ind==(max(Pxx_ind(fq,1))));
            [b,a]=butter(2,[(freq-5)/(0.5*samprate) (freq+5)/(0.5*samprate)],'bandpass');
%             thal_contact=[thal_contact r];
            bua_coh1=[bua_coh1; data{r,1}(rr,:)];
            ctx_coh(r,:)= data{r,1}(1,:);
            coherence1=[coherence1;Pxx_ind'];
            filt_thal1=[filt_thal1 ; filtfilt(b,a,data{r,1}(rr,:))]; %filt coherent subcortical contacts
            env_thal1=[env_thal1 ; abs(hilbert(filtfilt(b,a,data{r,1}(rr,:))))]; %filt coherent subcortical contacts
            phase_thal1=[phase_thal1 ; angle(hilbert(filtfilt(b,a,data{r,1}(rr,:))))]; %filt coherent subcortical contacts
            power_coh1=[power_coh1 ; power'];
            power_ctx(r,:)=power1;
           
            channel_bursts=data{r,1}(1,:);
            filt_ctx(r,:)=filtfilt(b,a,channel_bursts);
            env_ctx(r,:)=abs(hilbert(filtfilt(b,a,channel_bursts)));
            phase_ctx(r,:)=angle(hilbert(filtfilt(b,a,channel_bursts)));
            
    
        end
    end
    filt_thal{r,1}=filt_thal1;
    power_thal{r,1}=power_coh1;
    env_thal{r,1}=env_thal1;
    phase_thal{r,1}=phase_thal1;
    bua_coh{r,1}=bua_coh1;
    ctx_sub_coh{r,1}=coherence1;
end


l=[];
for i=1:size(phase_thal,1)
if ~isempty(phase_thal{i,1})
l=[l i];
end
end
phase_ctx = phase_ctx(any(phase_ctx,2),:);
coherence = coherence(any(coherence,2),:);
env_ctx = env_ctx(any(env_ctx,2),:);
filt_ctx = filt_ctx(any(filt_ctx,2),:);
power_ctx = power_ctx(any(power_ctx,2),:);
ctx_coh = ctx_coh(any(ctx_coh,2),:);

power_thal=power_thal(~cellfun('isempty',power_thal));
phase_thal=phase_thal(~cellfun('isempty',phase_thal));
env_thal=env_thal(~cellfun('isempty',env_thal));
filt_thal=filt_thal(~cellfun('isempty',filt_thal));

bua_coh=bua_coh(~cellfun('isempty',bua_coh));
ctx_sub_coh=ctx_sub_coh(~cellfun('isempty',ctx_sub_coh));
clearvars -except coherence env_ctx env_thal filt_ctx filt_thal  phase_ctx phase_thal power_thal power_ctx thal_contact samprate


CZ.animals= (A(lesion(n),1:33));
CZ.bua=data;
CZ.coh_animals=(A(lesion(n(l)),1:33));
CZ.env_thal=env_thal;
CZ.power_thal=power_thal;
CZ.phase_thal=phase_thal;
CZ.filt_thal=filt_thal;
CZ.env_ctx=env_ctx;
CZ.power_ctx=power_ctx;
CZ.phase_ctx=phase_ctx;
CZ.filt_ctx=filt_ctx;
CZ.samprate=samprate;
CZ.data_coh=data(l);

CZ.ctx_coh=ctx_coh;
CZ.bua_coh=bua_coh;
CZ.ctx_sub_coh=ctx_sub_coh;


clear all
% cd('C:\Users\creis\Documents\GitHub\CRcode\codes_thal')
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal')
load('SNR.mat')

for ik=1:size(SNR.env_ctx,1);
 env=SNR.env_ctx(ik,:);
 SNR.offset_raw{ik,1}=bursts(env);
 Ecogfiltered=SNR.filt_ctx(ik,:);
 SNR.offset_phase_al{ik,1}=bursts_aligned(env,Ecogfiltered);
end
    

for i=1:11
SNR.offset_raw_all{i,1}=sort(cell2mat(SNR.offset_raw{1,i}'),'ascend')
end


CZ.idrat=([1 7 9 11])';
SNR.idrat=([1 3 11 12])';
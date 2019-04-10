clear all
cd ('\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
load ('animal_lesion_nolesion.mat','A','lesion')
load 'SNr_ctx_probe.mat'


n=[];
for i=1:29
if size(WaveData_DCall{i,1},1)>1
n=[n i];
end
end
l=[];
for i=1:14
if ~isempty(phase_thal{i,1})
l=[l i];
end
end

clear all
close all
%  cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
load ('data_all')
load 'snr_ctx_probe'
fq=15:30;

samprate=1000;

[b,a]=butter(2,[15/(0.5*samprate) 30/(0.5*samprate)],'bandpass');
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
    m=1;
    [power1,f1]=pwelch(data{r,1}(1,:),1000,[],1000,1000);
    for rr=2:size(data{r,1},1) %contacts
        [power,f]=pwelch(data{r,1}(rr,:),1000,[],1000,1000); %ctx power spectracontact power spectra
        [Pxx_ind,F_ind]=mscohere(data{r,1}(1,:),data{r,1}(rr,:),samprate,[],samprate,samprate); %Magnitude-squared coherence between ctx-a given contact
        if  sum(Pxx_ind(fq))./sum(Pxx_ind(1:end))>0.1
%             thal_contact=[thal_contact r];
            coherence(r,:)=Pxx_ind;
            filt_thal1=[filt_thal1 ; filtfilt(b,a,data{r,1}(rr,:))]; %filt coherent subcortical contacts
            env_thal1=[env_thal1 ; abs(hilbert(filtfilt(b,a,data{r,1}(rr,:))))]; %filt coherent subcortical contacts
            phase_thal1=[phase_thal1 ; angle(hilbert(filtfilt(b,a,data{r,1}(rr,:))))]; %filt coherent subcortical contacts
            power_coh1=[power_coh1 ; power];
            power_ctx(r,:)=power1;
           
            channel_bursts=data{r,1}(1,:);
            filt_bursts= find(power1==(max(power1(fq,1))));
            [bb,aa]=butter(2,[(filt_bursts-5)/(0.5*samprate) (filt_bursts+5)/(0.5*samprate)],'bandpass');
            filt_ctx(r,:)=filtfilt(bb,aa,channel_bursts);
            env_ctx(r,:)=abs(hilbert(filtfilt(bb,aa,channel_bursts)));
            phase_ctx(r,:)=angle(hilbert(filtfilt(bb,aa,channel_bursts)));
            
    
        end
    end
    filt_thal{r,1}=filt_thal1;
    power_thal{r,1}=power_coh1;
    env_thal{r,1}=env_thal1;
    phase_thal{r,1}=phase_thal1;
end
for i=1:11
thal_coh(i,:)=data(thal_contact(i),:);
end

phase_ctx = phase_ctx(any(phase_ctx,2),:);
power_thal=power_thal(~cellfun('isempty',power_thal));

clearvars -except coherence env_ctx env_thal filt_ctx filt_thal  phase_ctx phase_thal power_thal power_ctx thal_contact samprate


SNR.animals= (A(lesion(n),1:33));
SNR.bua=data;
SNR.coh_animals=(A(lesion(n(l)),1:33));
SNR.env_thal=env_thal;
SNR.power_thal=power_thal;
SNR.phase_thal=phase_thal;
SNR.filt_thal=filt_thal;
SNR.env_ctx=env_ctx;
SNR.power_ctx=power_ctx;
SNR.phase_ctx=phase_ctx;
SNR.filt_ctx=filt_ctx;
SNR.samprate=samprate;
SNR.data_coh=data(l);


clear all
cd('C:\Users\creis\Documents\GitHub\CRcode\codes_thal')
cd('/Users/Carolina/Documents/GitHub/CRcode/codes_thal')
load('SNr.mat')

for ik=1:size(SNR.env_ctx,1);
 env=SNR.env_ctx(ik,:);
 SNR.ctxb_raw{1,ik}=bursts(env);
 Ecogfiltered=SNR.filt_ctx(ik,:);
 SNR.ctxb_phase_al{1,ik}=bursts_aligned(env,Ecogfiltered);
end
    

for i=1:14
SNR.ctxb_all{i,1}=sort(cell2mat(SNR.ctxb_raw{1,i}'),'ascend')
end

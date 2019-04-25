clear all
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal')
load('GPE data.mat')


n=[];
for i=1:29
    if size(WaveData_DCall{i,1},1)>1
        n=[n i];
    end
end

for i=1:length(n)
    r(i,:)=size(WaveData_DCall{n(i),1},2)
end
stp=find(r>=98000);
idx=n(stp);
for i=1:length(idx)
    data{i,1}=WaveData_DCall{idx(i),1}(:,1:98000);
end

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
        if  sum(Pxx_ind(15:35))./sum(Pxx_ind(1:end))>0.1
            freq= find(Pxx_ind==(max(Pxx_ind(15:35,1))));
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
BZ.phase_ctx = phase_ctx(any(phase_ctx,2),:);
BZ.env_ctx = env_ctx(any(env_ctx,2),:);
BZ.filt_ctx = filt_ctx(any(filt_ctx,2),:);
BZ.power_ctx = power_ctx(any(power_ctx,2),:);
BZ.ctx_coh = ctx_coh(any(ctx_coh,2),:);

BZ.power_thal=power_thal(~cellfun('isempty',power_thal));
BZ.phase_thal=phase_thal(~cellfun('isempty',phase_thal));
BZ.env_thal=env_thal(~cellfun('isempty',env_thal));
BZ.filt_thal=filt_thal(~cellfun('isempty',filt_thal));

BZ.bua_coh=bua_coh(~cellfun('isempty',bua_coh));
BZ.ctx_sub_coh=ctx_sub_coh(~cellfun('isempty',ctx_sub_coh));


clearvars -except BZ

for ik=1:size(BZ.env_ctx,1);
    env=BZ.env_ctx(ik,:);
    BZ.onset_raw{1,ik}=bursts(env);
    BZ.offset_raw{1,ik}=bursts_off(env);
    Ecogfiltered=BZ.filt_ctx(ik,:);
    BZ.onset_phase_al{1,ik}=bursts_aligned(env,Ecogfiltered);
    BZ.offset_phase_al{1,ik}=bursts_aligned_off(env,Ecogfiltered);
end


for i=1:size(BZ.env_ctx,1)
    BZ.onset_raw_all{i,1}=sort(cell2mat(BZ.onset_raw{1,i}'),'ascend')
    BZ.offset_raw_all{i,1}=sort(cell2mat(BZ.offset_raw{1,i}'),'ascend')
    BZ.off_phaseal_all{i,1}=sort(cell2mat(BZ.offset_phase_al{1,i}'),'ascend')
    BZ.on_phaseal_all{i,1}=sort(cell2mat(BZ.onset_phase_al{1,i}'),'ascend')
end


% % clearvars -except BZ
% % cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal')
% % save 'BZ'





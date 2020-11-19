%a)select the region we want to analyse 
%b)do the spike sorting - BUA
%c)filter the EGG in the beta range
%d)leave thal region unfiltered

clear all
cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A2_Thal/MAT')
%  cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A2_Thal\mat')

load 'animal_lesion_nolesion.mat'

state_mat=lesion;% lesion vs nolesioned
nfile_i=3:(newfile);
rr=1;
for nn=1:length(state_mat);
    
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/MATLAB/KN')
%    cd ('\Users\creis\Documents\MATLAB\KN')
    clearvars -except nfile_i state_mat A rr  WaveDatall WaveData_DCall evoked_all ctx_all nn coh_all beta_coh_all power_all onset_all
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
    
    if (length ([thal_BZ thal_CZ thal_ZI thal_AV thal_Rt thal_AM thal_SNr thal_NW thal_PC]))== length (thalchan)
        thal_local=thalchan(thal_CZ); %%% Chose thal location
    else thal_local=newchannel;
    end
    
%     nu=[];
%     mua=[];% also has mua/bua cases on them
%     sua=[];
%     bua=[];
%     for i=1:length(thal_local);
%         eval(['note=' B{thal_local(i)} '.note_labeledcell']);
%         if ~isempty(min(find(note=='n')) & min(find(note=='o')) & min(find(note=='u')) & min(find(note=='s')) & min(find(note=='e')));
%             nu=[nu,i];
%         elseif ~isempty(min(find(note=='M')) & min(find(note=='U')) & min(find(note=='A')));
%             mua=[mua,i];
%         elseif ~isempty(min(find(note=='S')) & min(find(note=='U')) & min(find(note=='A')));
%             sua=[sua,i];
%         elseif ~isempty(min(find(note=='B')) & min(find(note=='U')) & min(find(note=='A')));
%             bua=[bua,i];
%         end
%     end
%     chanofinterest= thal_local([mua sua bua nu]);
%     %   length(thal_local)==length([mua bua sua nu]) % check all channels included
    
    chanofinterest= thal_local;
    eval(['samprateold=1/' B{ctxchan} '.interval;']);
    eval(['WaveData(1,:)=' B{ctxchan} '.values;']);
    WaveData=double(WaveData);
   
    cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A2_Thal/code')

%     cd ('\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A2_Thal\code')

    ts=timeseries(WaveData(1,:),0:(1/samprateold):((size(WaveData,2)-1)/samprateold));
    ts1=resample(ts,0:0.001:((size(WaveData,2)-1)/samprateold),'linear');
    WaveData_DC(1,:)=ts1.data;
    timeold=0:(1/samprateold):(length(WaveData(1,:))-1)*(1/samprateold);
    for ii=1:length(chanofinterest)
        eval(['WaveData(1+ii,:)=' B{chanofinterest(ii)} '.values;']);
        WaveData=double(WaveData);
        WaveData_DC(1+ii,:)=makemua_CR_1(WaveData(1+ii,:)',0.001,0.003,round(samprateold),1000,4);
    end
    thal_raw= WaveData_DC(2:(end),:);
    
    clear muafilt
    samprate=1000;
    time=0:(1/samprate):(length(WaveData_DC(1,:))-1)*(1/samprate);
    muafilt=[];
    [Pxx_ind,F_ind]=pwelch(WaveData_DC(1,:),samprate,[],samprate,samprate);
    Pxx_ind_beta=Pxx_ind(16:36);
    filtrange=14+find(Pxx_ind_beta==max(Pxx_ind_beta));
    [b,a]=butter(2,[(filtrange-5)/(0.5*samprate) (filtrange+5)/(0.5*samprate)],'bandpass');
    muafilt(1,1:length(time)) = filtfilt(b,a,WaveData_DC(1,:));
    
    run 'phase_onset_new.m'
    coh_one=[];
    beta_coh=[];
    power_probe=[];
    ind_p=[];
    m=2;
    n=1;
    for r=2:size(WaveData_DC,1)
        [power,f]=pwelch(WaveData_DC(r,:),1000,[],1000,1000);
        [Pxx_ind,F_ind]=mscohere(WaveData_DC(1,:),WaveData_DC(r,:),samprate,[],samprate,samprate);
        Pxx_ind_beta=Pxx_ind(16:36);
        filtrange=14+find(Pxx_ind_beta==max(Pxx_ind_beta));
        [b,a]=butter(2,[(filtrange-5)/(0.5*samprate) (filtrange+5)/(0.5*samprate)],'bandpass');
        P=sum(Pxx_ind(16:36))./sum(Pxx_ind(1:end));
        power_probe(n,:)=(power)';
        coh_one(n,:)=(Pxx_ind)';
        beta_coh(:,n)=P;
        n=n+1;
    end
    
    power_all{rr,:}=power_probe;
    coh_all{rr,:}=coh_one;
    beta_coh_all{rr,:}=beta_coh;
    if ~isnan (out_evoked)
        evoked_all{rr,1}= out_evoked;
        ctx_all{rr,1}=output_ctx;
        %onset_all(rr,:)=[size(output,1)*size(output,2)]; Uncomment and sum
        %to get nu,mber of trials (contacts x  onset points)
    end
     onset_all(rr,:)=length(onset);
     WaveData_DCall{rr,:}= WaveData_DC;
WaveDatall{rr,:}=WaveData;
    rr=rr+1;
   
end

% evoked=evoked_all(~cellfun('isempty',evoked_all));
% clearvars -except power_all evoked_all power_all beta_coh_all coh_all evoked
% cd('\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A2_Thal\mat')
% save 'lesion_param_2000_SNr_new.mat'
% run 'new_thal.m'

WDC=WaveData_DCall
for j=1:size(WDC,1)
    if  size(WDC{j},1)==1
        WDC{j}=[];
    end
end

 data= WDC(~cellfun('isempty', WDC));
 recs=A(state_mat([find(~cellfun('isempty',WDC))]),:);


clearvars -except WaveData_DCall data
% cd('\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
% save 'CZ_ctx_probe.mat'

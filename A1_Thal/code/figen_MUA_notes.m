
clear all
%cd ('C:\Users\wolf5173\Documents\Carolina_code\codes_thal\A1_Thal\mat')
cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A1_Thal/mat')
load 'animal_lesion_nolesion.mat'
newfile=size(A,1);
state_mat=lesion;% lesion vs nolesioned
nfile_i=3:(newfile);
med_thal_cell=cell(1,1);
med_ctx=[];
stats_bursts=cell(1,1);
x=1;
check=[];
for nn=3;
    %1:length(state_mat);
    
    %cd C:\Users\wolf5173\Documents\MATLAB\KN
    cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/MATLAB/KN')
    clearvars -except x med_thal_cell med_ctx nfile_i state_mat nn A check stats_bursts
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
        thal_local=thalchan(thal_SNr); %%% Chose thal location
    else thal_local=newchannel;
    end
    
    
    nu=[];
    mua=[];% also has mua/bua cases on them
    sua=[];
    bua=[];
    for i=1:length(thal_local);
        eval(['note=' B{thal_local(i)} '.note_labeledcell']);
        if ~isempty(min(find(note=='n')) & min(find(note=='o')) & min(find(note=='u')) & min(find(note=='s')) & min(find(note=='e')));
            nu=[nu,i];
        elseif ~isempty(min(find(note=='M')) & min(find(note=='U')) & min(find(note=='A')));
            mua=[mua,i];
        elseif ~isempty(min(find(note=='S')) & min(find(note=='U')) & min(find(note=='A')));
            sua=[sua,i];
        elseif ~isempty(min(find(note=='B')) & min(find(note=='U')) & min(find(note=='A')));
            bua=[bua,i];
        end
    end
    
    
    chanofinterest= thal_local([mua sua bua nu]);
    
%   length(thal_local)==length([mua bua sua nu]) % check all channels included
    
    
    eval(['samprateold=1/' B{ctxchan} '.interval;']);
    eval(['WaveData(1,:)=' B{ctxchan} '.values;']);
    WaveData=double(WaveData);
    
%     cd C:\Users\wolf5173\Documents\Carolina_code\codes_thal\A1_Thal\code
cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A1_Thal/code')
    
    ts=timeseries(WaveData(1,:),0:(1/samprateold):((size(WaveData,2)-1)/samprateold));
    ts1=resample(ts,0:0.001:((size(WaveData,2)-1)/samprateold),'linear');
    WaveData_DC(1,:)=ts1.data;
    timeold=0:(1/samprateold):(length(WaveData(1,:))-1)*(1/samprateold);
    for ii=1:length(chanofinterest)
        eval(['WaveData(1+ii,:)=' B{chanofinterest(ii)} '.values;']);
        WaveData=double(WaveData);
        WaveData_DC(1+ii,:)=makemua_CR_1(WaveData(1+ii,:)',0.001,0.003,round(samprateold),1000,4);
    end
    
    clear muafilt
    samprate=1000;
    
    time=0:(1/samprate):(length(WaveData_DC(1,:))-1)*(1/samprate);
    
    muafilt=[];
    [Pxx_ind,F_ind]=pwelch(WaveData_DC(1,:),samprate,[],samprate,samprate);
    Pxx_ind_beta=Pxx_ind(16:36);
    filtrange=14+find(Pxx_ind_beta==max(Pxx_ind_beta));
    [b,a]=butter(2,[(filtrange-5)/(0.5*samprate) (filtrange+5)/(0.5*samprate)],'bandpass');
    muafilt(1,1:length(time)) = filtfilt(b,a,WaveData_DC(1,:));
    
    ind_p=[];
    m=2;
    for r=2:size(WaveData_DC,1)
        [Pxx_ind,F_ind]=mscohere(WaveData_DC(1,:),WaveData_DC(r,:),samprate,[],samprate,samprate);
        Pxx_ind_beta=Pxx_ind(16:36);
        filtrange=14+find(Pxx_ind_beta==max(Pxx_ind_beta));
        [b,a]=butter(2,[(filtrange-5)/(0.5*samprate) (filtrange+5)/(0.5*samprate)],'bandpass');
        P=sum(Pxx_ind(16:36))./sum(Pxx_ind(1:end));
        if P>0.1
            muafilt(m,1:length(time)) = filtfilt(b,a,WaveData_DC(r,:));
            m=m+1;
            ind_p=[ind_p,r];
        end
    end
    
    
    %     CHECK1
    
    run 'median_bursts.m'
    
    
    if ~isnan (med_thal_sub)
        med_thal_cell{x,1}= med_thal_sub;
        stats_bursts{x,1}=[median_b SD_b];
        x=x+1;
    end
    med_ctx=[med_ctx;med_ctx_sub];
    med_ctx( ~any(med_ctx,2), : ) = [];
    
%     figure(1)
%     for  i = (1:size(WaveData_DC,1))
%         [pxx,f]=pwelch(WaveData_DC(i,:),1000,[],1000,1000);
%         plot(pxx(5:100))
%         hold on
%     end
%     
%     if ~isempty (ind_p)
%         for f=ind_p(1):ind_p (end)
%             [pxx,f]=pwelch(WaveData_DC(f,:),1000,[],1000,1000);
%             figure(2)
%             plot(pxx(5:100),'r.')
%         end
%     end  
end

med_thal=cell2mat([med_thal_cell]);

clearvars -except med_thal med_ctx WaveData_DC time stats_bursts

% 
figure(2)
subplot (2,1,1)
time_plot=[-500:1:500];
plot(time_plot,median(med_ctx))
hold on
plot(time_plot,prctile(med_ctx,75));
plot(time_plot,prctile(med_ctx,25));
subplot(2,1,2)
time_plot=[-500:1:500];
plot(time_plot,median(med_thal))
hold on
plot(time_plot,prctile(med_thal,75));
plot(time_plot,prctile(med_thal,25));


clear all
% cd ('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\probe SUA_act_mat')
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/probe SUA_act_mat')
load list_probe_sua.mat; A=probe_list; clearvars probe_list;
newfile=[];
nolesion=[];
lesion=[];
newfile=size(A,1);
for nfile_i=1:(newfile)
    name=A(nfile_i,1:(find(A(nfile_i,:)=='.')-1));
    load(name)
    if IpsiEEG.dopamine=='6-OHDA'
        lesion=[lesion nfile_i];
        
    else nolesion=[nolesion nfile_i];
        
    end
    clearvars -except A nfile new file lesion nolesion
end

state_mat=lesion;% lesion vs nolesioned
newfile=size(A,1);

for nn=1:length(state_mat);  %%% chose state (lesion OR nolesion)
    clearvars -except nfile_i state_mat A data_all Ecog_all nn
    
    name=A(state_mat(nn),1:(find(A((state_mat((nn))),:)=='.')-1))
    load(name)
    B=who;
    
    ctxchan=[];
    thalchan=[];
    for ii=1:size(B,1)
        if ~isempty(min(find(B{ii}=='I')) & min(find(B{ii}=='p')) & min(find(B{ii}=='i')) &  min(find(B{ii}=='E')) & min(find(B{ii}=='E')));
            ctxchan=ii;
        elseif ~isempty(min(find(B{ii}=='p')) & min(find(B{ii}=='r')) & min(find(B{ii}=='o')) & min(find(B{ii}=='b')));
            thalchan=[thalchan ii];
        end
    end
    
    thal_BZ=[];%basal-ganglia recipient thalamus
    thal_PF=[];
    thal_ZI=[];%Zona incerta
    thal_AV=[];%Amteroventral thalamus
    thal_SNr=[];%Substantia Nigra reticulata
    thal_NW=[];%No where
    thal_PC=[];%Para-central thalamus
    thal_VPPC=[];%Antero medial thalamus
    thal_CZ=[];%cerebellar recipient thalamus
    thal_CL=[];
    thal_VPM=[];
    for i=1:length(thalchan)
        eval(['location=' B{thalchan(i)} '.location']);
        if ~isempty(min(find(location=='B')) & min(find(location=='Z')));
            thal_BZ=[thal_BZ,i];
            %         elseif ~isempty(min(find(location=='P')) & min(find(location=='F')));
            %             thal_PF=[thal_PF,i];
            %         elseif ~isempty(min(find(location=='A')) & min(find(location=='V')));
            %             thal_AV=[thal_AV,i];
            %         elseif ~isempty(min(find(location=='Z')) & min(find(location=='I')));
            %             thal_ZI=[thal_ZI,i];
        elseif ~isempty(min(find(location=='S')) & min(find(location=='N')) & min(find(location=='r')));
            thal_SNr=[thal_SNr,i];
            %         elseif ~isempty(min(find(location=='n')) & min(find(location=='o')) & min(find(location=='w')));
            %             thal_NW=[thal_NW,i];
            %         elseif ~isempty(min(find(location=='P')) & min(find(location=='C')));
            %             thal_PC=[thal_PC,i];
            %         elseif ~isempty(min(find(location=='V')) & min(find(location=='P')) & min(find(location=='P')) & min(find(location=='C')));
            %             thal_VPPC=[thal_VPPC,i];
        elseif ~isempty(min(find(location=='C')) & min(find(location=='Z')));
            thal_CZ=[thal_CZ,i];
            %         elseif ~isempty(min(find(location=='C')) & min(find(location=='L')));
            %             thal_CL=[thal_CL,i];
            %         elseif ~isempty(min(find(location=='V')) & min(find(location=='P'))& min(find(location=='M')));
            %             thal_VPM=[thal_VPM,i];
        end
    end
    
    %     if (length ([thal_BZ thal_VPM thal_PF thal_ZI thal_AV thal_PC thal_CZ thal_VPPC thal_SNr thal_NW thal_PC]))== length (thalchan)
    thal_local=thalchan(thal_CZ); %%% Chose thal location
    %     else
    %         thal_local=[];
    %     end
    
    eval(['samprateold=1/' B{ctxchan} '.interval;']);
    eval(['WaveData=' B{ctxchan} '.values;']);
    WaveData=double(WaveData)';
    ts=timeseries(WaveData(1,:),0:(1/samprateold):((size(WaveData,2)-1)/samprateold));
    ts1=resample(ts,0:0.001:((size(WaveData,2)-1)/samprateold),'linear');
    WaveData_DC(1,:)=ts1.data;
    
    for il=1:length(thal_local);
        eval(['sr=1/' B{thal_local(il)} '.interval;']);
        
        srn=1000;
        eval(['dataold=' B{thal_local(il)} '.values;']);
        dataold=full(dataold);
        data=zeros(1,100000);
        timeold=0:1/sr:(size(dataold,1)-1)/sr;
        time=0:1/srn:(size(data,2)-1)/srn;
        spk_t=timeold(find(dataold==1));
        spk_tround=round(spk_t,3);
        nnn=[];
        for ff=1:length(spk_t)
            [ d, ix ] = min( abs( time-spk_tround(ff) ) );
            nnn=[nnn ix];
        end
        data(nnn)=1;
        data_all{nn,1}(il,:)=data;
        Ecog_all(nn,:)=WaveData_DC;
        %         data_ones=find(data==1);
    end
end

data_region=data_all(~cellfun('isempty',data_all))
Ecog_region=Ecog_all(any(Ecog_all,2),:);


cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/probe SUA_act_mat')
clearvars -except data_all Ecog_all A lesion time data_region Ecog_region
save 'data_SUA_CZ'







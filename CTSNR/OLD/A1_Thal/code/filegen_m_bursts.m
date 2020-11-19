clear all
cd C:\Users\wolf5173\Documents\MATLAB\KN
A=ls;
newfile=size(A,1);
cd ('C:\Users\wolf5173\Documents\Carolina_code\codes_thal\A1_Thal\mat')
load 'animal_lesion_nolesion.mat'
load 'thal_BZ_all.mat'

state_mat=lesion;% lesion vs nolesioned
nfile_i=3:(newfile);
med_thal_cell=cell(1,1);
med_ctx=[];
x=1;
check=[];
thal_all=cell(1,1);
for nn=1:length(subjects);
    
    cd C:\Users\wolf5173\Documents\MATLAB\KN
    clearvars -except x med_thal_cell med_ctx nfile_i state_mat nn A check subjects
    name=A(state_mat(nn),1:(find(A((state_mat((nn))),:)=='.')-1))
    load(name)
    B=who;
    
    ctxchan=[];
    thalchan=[];
    for ii=1:size(B,1)
        if ~isempty(min(find(B{ii}=='i')) & min(find(B{ii}=='p')) & min(find(B{ii}=='s')))
            ctxchan=ii;
        end
    end
    
    cd ('C:\Users\wolf5173\Documents\Carolina_code\codes_thal\A1_Thal\mat')
    load 'thal_BZ_all.mat'
    chanofinterest= thal_all{nn};
    
%   length(thal_local)==length([mua bua sua nu]) % check all channels included
    
    
    eval(['samprateold=1/' B{ctxchan} '.interval;']);
    eval(['WaveData(1,:)=' B{ctxchan} '.values;']);
    WaveData=double(WaveData);
    
    cd C:\Users\wolf5173\Documents\Carolina_code\codes_thal\A1_Thal\code
    
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

    run 'median_bursts.m'

    if ~isnan (med_thal_sub)
        med_thal_cell{x,1}= med_thal_sub;
        x=x+1;
    end
    med_ctx=[med_ctx;med_ctx_sub];
    med_ctx( ~any(med_ctx,2), : ) = [];
     
end

med_thal=cell2mat([med_thal_cell]);
% 
% clearvars -except med_thal med_ctx WaveData_DC time
% 
% % 
% figure(2)
% subplot (2,1,1)
% time_plot=[-500:1:500];
% plot(time_plot,median(med_ctx))
% hold on
% plot(time_plot,prctile(med_ctx,75));
% plot(time_plot,prctile(med_ctx,25));
% subplot(2,1,2)
% time_plot=[-500:1:500];
% plot(time_plot,median(med_thal))
% hold on
% plot(time_plot,prctile(med_thal,75));
% plot(time_plot,prctile(med_thal,25));
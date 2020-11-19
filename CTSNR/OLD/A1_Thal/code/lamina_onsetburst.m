clear all

cd C:\Users\wolf5173\Documents\Carolina_code\codes_thal\A2_Thal\mat
load 'animal_lesion_nolesion'

state_mat=lesion;% lesion vs nolesioned
nfile_i=3:(newfile);
med_thal_cell=cell(1,1);
med_ctx=[];
x=1;
n=1;
check=[];
for nn=1:length(state_mat);
    
    cd C:\Users\wolf5173\Documents\MATLAB\KN
    clearvars -except x med_thal_cell med_ctx nfile_i state_mat nn A n median_indv beta_coh_all
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
    
    thal_BZ=[];
    for i=1:length(thalchan)
        eval(['location=' B{thalchan(i)} '.location']);
        if ~isempty(min(find(location=='B')) & min(find(location=='Z')));
            thal_BZ=[thal_BZ,i];
        end
    end
    
    if ~isempty (thal_BZ)
        chanofinterest=thalchan(thal_BZ);
        
        eval(['samprateold=1/' B{ctxchan} '.interval;']);
        eval(['WaveData(1,:)=' B{ctxchan} '.values;']);
        WaveData=double(WaveData);
        
        cd C:\Users\wolf5173\Documents\Carolina_code\codes_thal\A2_Thal\code
        
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
            muafilt(m,1:length(time)) = filtfilt(b,a,WaveData_DC(r,:));
            beta_coh(:,m)=P;
            m=m+1;
        end
        
        cd C:\Users\wolf5173\Documents\Carolina_code\codes_thal
        run 'processing_individuals.m'
        
        median_indv{n,1}=median_chan;
        beta_coh_all{n,1}=beta_coh;
        n=n+1;
       
    end
end



for f=1:size(median_indv,1)
close all
power_probe=median_indv{f,1};
chan=1:size(power_probe,1);
time=[-500:1:500];
d1 = 1:length(chan);
d2 = 1:length(time)
Hc_all(d1,d2) = power_probe;  

% freq_band=[15:35];
% beta_all=[];
% for b=1:length(chan)
% beta_all(b,1:length(freq_band)) =Hc_all(b,(freq_band));
% norm_all(b)=sum(Hc_all(b,1:501));
% norm_beta(b,1:length(freq_band))=(100*beta_all(b,:)-norm_all(b))./norm_all(b);
% end
% 
% plot_beta=norm_beta;


figure(1)
imagesc(time,chan,Hc_all)
axis square; colorbar
xlim([-500 500]); 
ylim([0 size(power_probe,1)+1]);

% figure(2)
% for c=1:size(power_probe,1)
% subplot(size(power_probe,1),1,c)
% plot(time,power_probe(c,:))
% hold on
% end

end



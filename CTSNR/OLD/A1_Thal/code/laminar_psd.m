clear all

cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A1_Thal/mat')
load 'animal_lesion_nolesion'

state_mat=lesion;% lesion vs nolesioned
nfile_i=3:(newfile);
beta_coh_all=cell(1,1);
coh_all=cell(1,1);
power_all=cell(1,1);
x=1;
check=[];
rr=1;
for nn=1:length(state_mat)
    
    cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/MATLAB/KN')
    clearvars -except x med_thal_cell med_ctx nfile_i state_mat nn A check coh_all beta_coh_all rr power_all
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
    
        thal_SNR=[];
    for i=1:length(thalchan)
        eval(['location=' B{thalchan(i)} '.location']);
        if ~isempty(min(find(location=='S')) & min(find(location=='N')) &  min(find(location=='r')));
            thal_SNR=[thal_SNR,i];
        end
    end
    
    
    thal_local=thal_SNR;

    if ~isempty (thal_local)
        chanofinterest=thalchan(thal_local);
        
        eval(['samprateold=1/' B{ctxchan} '.interval;']);
        eval(['WaveData(1,:)=' B{ctxchan} '.values;']);
        WaveData=double(WaveData);
        
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
            
%             max(Pxx_ind(10:50))
%             F_ind(find(Pxx_ind==max(Pxx_ind(10:50))))
            
            power_probe(n,:)=(power)';
            coh_one(n,:)=(Pxx_ind)';
            beta_coh(:,n)=P;
            n=n+1;
        end
        
        power_all{rr,:}=power_probe;
        coh_all{rr,:}=coh_one;
        beta_coh_all{rr,:}=beta_coh;
        rr=rr+1;
    end
    
end

clearvars -except power_all coh_all beta_coh_all



% for i=1:size(power_all,1)
% power_probe=power_all{i,1};
% 
% chan=1:size(power_probe,1);
% freq= (f)';
% d1 = 1:length(chan);
% d2 = 1:length(freq);
% Hc_all(d1,d2) = power_probe;  % check check and check this !
% 
% 
% freq_band=[5:50];
% beta_all=[];
% for b=1:length(chan)
% beta_all(b,1:length(freq_band)) =Hc_all(b,(freq_band));
% norm_all(b)=sum(Hc_all(b,1:501));
% norm_beta(b,1:length(freq_band))=(100*beta_all(b,:)-norm_all(b))./norm_all(b);
% end
% 
% plot_beta=norm_beta;
% 
% figure (1)
% imagesc(freq(freq_band)+1,chan,plot_beta)
% axis square; colorbar
% xlim([5 50]); 
% ylim([1 size(power_probe,1)]);
% 
% figure(2)
% for c=1:size(power_probe,1)
% % subplot(size(power_probe,1),1,c)
% plot(freq_band,power_probe(c,freq_band))
% hold on
% legend('show')
% end
% 
% end
% 
% 
% 
% 

% clearvars -except beta_coh_all coh_all power_all WaveData_DC 
% 
% for i=1:length(beta_coh_all)
%     plot(beta_coh_all{i},'o')
%     hold on
%     plot([0,12],[0.1,0.1])
% end
% 
% 
% beta_all=power_all(:,15:36);
% for i=1:size(beta_all,1)
%     max_bchan(i,1)= (max(beta_all(i,:)));
% end
% 
% for i=1:size(beta_all,1)
%     freq_max(i,1)=find(beta_all(i,:)==max_bchan(i));
% end


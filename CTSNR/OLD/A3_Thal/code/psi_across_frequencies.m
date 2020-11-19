clear all
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
% cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
A=char('AV_ctx_probe.mat','BZ_ctx_probe.mat ','CZ_ctx_probe.mat ','SNR_ctx_probe.mat','Rt_ctx_probe.mat ');
samprate=1000;
% % load 'bz_ctx_probe'; color_b=[0.5 0 0.5];
% load 'snr_ctx_probe'; color_b=[0 0 0.5];

colors = { [0.2 0.3 0.4]; [0.5 0 0.5] ; [0.5 0.5 0.5]; [0 0 0.5]; [0.8 0.8 0.8]};
%
freq=[6:3:48];
% freq=[15:5:25];

for xx=1:size(A,1);
    name=A(xx,1:(find(A(xx,:)=='.')-1));
    load(name)
    color_b=colors{xx,1};
    clearvars -except xx colors A freq samprate data color_b psi_all
    for t=1:size(freq,2)
        psi_rec=[];
        ctcts=0;
        for r=1:size(data,1); %electrodes
            coh_thal=[];
            clearvars a b euler
            pp=1;
            for rr=2:size(data{r,1},1) %contacts
                [power,f]=pwelch(data{r,1}(rr,:),1000,[],1000,1000); %ctx power spectracontact power spectra
                [b,a]=butter(2,[(freq(t)-1)/(0.5*samprate) (freq(t)+1)/(0.5*samprate)],'bandpass');
                [Pxx_ind,F_ind]=mscohere(data{r,1}(1,:),data{r,1}(rr,:),samprate,[],samprate,samprate); %Magnitude-squared coherence between ctx-a given contact
                euler1=[];
                %              if   sum(Pxx_ind(freq(t)-1:freq(t)+1))/sum(Pxx_ind(1:end))>0.1; %sum(power(freq(t)-5:freq(t)+5))./sum(power(1:end))>0.1
                ctcts=ctcts+1;
                coh_thal=[coh_thal sum(Pxx_ind(freq(t)-1:freq(t)+1))./sum(Pxx_ind(1:end))];
                non_norm=squeeze(angle(hilbert(filtfilt(b,a,data{r,1}(1,:))))-angle(hilbert(filtfilt(b,a,data{r,1}(rr,:))))); %dif angles between ctx-subctx
                
                for x =1:size(non_norm,2)
                    if non_norm(1,x)>pi
                        non_norm(1,x)=(non_norm(1,x))-(2.*pi);
                    elseif non_norm(1,x)<-pi
                        non_norm(1,x)=(non_norm(1,x))+(2.*pi);
                    else
                        non_norm(1,x)= non_norm(1,x);
                    end
                end
                
                dif_angs=non_norm;
                clearvars non_norm;
                euler1=abs(sum(exp(sqrt(-1)*(dif_angs)))./(length(dif_angs)));
                euler(pp,:)=euler1;
                pp=pp+1;
                
                psi_rec(r,:)=mean(euler);
                clearvars dif_angs
                %             end
            end
            
        end
        
        if ~isempty (psi_rec)
            psi_freq(t,:)=psi_rec;
            psi_all(xx,t)=mean(psi_rec);
        else
            psi_all(xx,t)=NaN;
        end
        ctcts_all(t,1)=ctcts;
    end
    
end


n=[];
for dum=1;
    for xx=1:size(A,1);
        color_b=colors{xx,1};
        p1=set(plot(freq,psi_all(xx,:)),'LineWidth',1.5,'Color',color_b);
        legend([p1],{'AV-CTX','BZ-CTX','CZ-CTX','SNr-CTX','Rt-CTX'})
        hold on
        n(1,xx)=mean(psi_all(xx,:))+std(psi_all(xx,:))
    end
%     th=mean(n);
%     thr=zeros(1,15);
%     thr(1,:)=th;
%     p1=set(plot(freq,thr,'--'),'LineWidth',1.5,'Color','black');
%     legend([p1],{'avg(mean+1std)'})
end



% xticklabels({'5-7','8-10','11-13','14-16','18-19','20-22','24-25','26-28','29-31','32-34','35-37','38-40','41-43','44-46','47-49'})
% legend ('VA-CTX','BZ-CTX','CZ-CTX','SNR-CTX','RT-CTX',"FontSize", 12)
xlabel ('Frequencies(Hz)',"FontSize", 12)
ylabel ('Phase Synchrony Index',"FontSize", 12)
box off
% for i=1:size(psi_rec,1)
% subplot(size(psi_rec,1),1,i)
% bar([psi_rec{i,1}])
% end

% figure()
% for ii=1:3
% for i=1:4
%     figure(ii)
%     subplot(1,size(ang_all,1),i)
%     polarhistogram(ang_all{i,ii},12)
% end
% end


% cd('/Users/Carolina/Desktop/TC_data')
% savefig('BUA_CTX_RT')
% saveas(gca,'BUA_CTX_RT.png')


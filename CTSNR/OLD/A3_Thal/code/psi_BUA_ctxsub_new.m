clear all
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
% cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
A=char('AV_ctx_probe.mat','BZ_ctx_probe.mat ','CZ_ctx_probe.mat ','SNR_ctx_probe.mat','Rt_ctx_probe.mat ');
samprate=1000;
colors = { [0.2 0.3 0.4]; [0.5 0 0.5] ; [0.5 0.5 0.5]; [0 0 0.5]; [0.8 0.8 0.8]};

freq={[8:13];[15:20];[21:30]} ;
for xx=1:size(A,1);
    name=A(xx,1:(find(A(xx,:)=='.')-1));
    load(name)
    color_b=colors{xx,1};
    clearvars -except xx colors A freq samprate data color_b psi_all ctcts_all
    for t=1:size(freq,1)
        psi_rec=[];
        ctcts=0;
        for r=1:size(data,1); %electrodes
            coh_thal=[];
            clearvars a b euler
            pp=1;
            for rr=2:size(data{r,1},1) %contacts
                [power,f]=pwelch(data{r,1}(rr,:),1000,[],1000,1000); %ctx power spectracontact power spectra
                [b,a]=butter(2,[(freq{t,1}(1))/(0.5*samprate) (freq{t,1}(end))/(0.5*samprate)],'bandpass');
%                 if t==length(freq)
%                     [b,a]=butter(2,[49/(0.5*samprate) 100/(0.5*samprate)],'bandpass');
%                 end
                [Pxx_ind,F_ind]=mscohere(data{r,1}(1,:),data{r,1}(rr,:),samprate,[],samprate,samprate); %Magnitude-squared coherence between ctx-a given contact
                euler1=[];
                if   sum(Pxx_ind(freq{t,1}))/sum(Pxx_ind(1:end))>0.1; %sum(power(freq(t)-5:freq(t)+5))./sum(power(1:end))>0.1
                    ctcts=ctcts+1;
                    coh_thal=[coh_thal sum(Pxx_ind(freq{t,1}))./sum(Pxx_ind(1:end))];
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
                end
            end
            
        end
        
        if ~isempty (psi_rec)
            psi_all(xx,t)=mean(psi_rec);
        else
            psi_all(xx,t)=NaN;
        end
        ctcts_all(xx,t)=ctcts;
    end
end

p1=bar(psi_all');
% p1.FaceColor=colors;
% xticklabels({'5-15','16-26','27-37','38-48','49-100'})
legend ('VA-CTX','BZ-CTX','CZ-CTX','SNR-CTX','RT-CTX',"FontSize", 12)
xticklabels({'{\alpha}(8-14)','low{\beta}(15-20)','high{\beta}(21-30)'})
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


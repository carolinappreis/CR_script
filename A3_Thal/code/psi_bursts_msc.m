clear all

cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
load ('data_all' , 'freq')
load 'BZ_ctx_probe'

samprate=1000;
time=0:(1/samprate):(length(data{1,1}(1,:))-1)*(1/samprate);

for t=1:size(freq,1)
    ctx_b1=[];
    clust_b=[];
    clust_s=[];
    all_psi_b=[];
    all_psi_sur=[];
    ctcts=0;
    for r=1:size(data,1); %electrodes
        thal_contact=[];
        coh_thal=[];
        filt_thal=[];
        clearvars a b
        %         plot(f,power1)
        %         xlim ([0 100])
        %         ylim ([0 200000])
        pp=1;
        for rr=2:size(data{r,1},1) %contacts
            [power,f]=pwelch(data{r,1}(rr,:),1000,[],1000,1000); %ctx power spectracontact power spectra
            %                         plot(f,power)
            %                         xlim ([0 100])
            
            
            [b,a]=butter(2,[(freq(t)-5)/(0.5*samprate) (freq(t)+5)/(0.5*samprate)],'bandpass');
            if t==length(freq)
                [b,a]=butter(2,[49/(0.5*samprate) 60/(0.5*samprate)],'bandpass');
            end
            [Pxx_ind,F_ind]=mscohere(data{r,1}(1,:),data{r,1}(rr,:),samprate,[],samprate,samprate); %Magnitude-squared coherence between ctx-a given contact
            if  sum(Pxx_ind(freq(t)-5:freq(t)+5))./sum(Pxx_ind(1:end))>0.1; %sum(power(freq(t)-5:freq(t)+5))./sum(power(1:end))>0.1
                ctcts=ctcts+1;
                thal_contact=[thal_contact rr];
                coh_thal=[coh_thal sum(Pxx_ind(freq(t)-5:freq(t)+5))./sum(Pxx_ind(1:end))];
                filt_thal=[filt_thal ; filtfilt(b,a,data{r,1}(rr,:))]; %filt coherent subcortical contacts
                
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
                dif_angs(rr,:)=non_norm;
                clearvars non_norm;
                
                    channel_bursts=data{r,1}(rr,:);
%                 channel_bursts=data{r,1}(1,:);
                
                [power1,f]=pwelch(channel_bursts,1000,[],1000,1000); %ctx power spectra
                filt_bursts= find(power1==(max(power1(13:30,1))));
                [b,a]=butter(2,[(filt_bursts-5)/(0.5*samprate) (filt_bursts+5)/(0.5*samprate)],'bandpass');
                env=abs(hilbert(filtfilt(b,a,channel_bursts)));
                threshold=prctile(env,75);
                tt(size(env,1):size(env,2))=threshold;
                indexexceed=find(env>threshold);
                diffindex=diff(indexexceed);
                pnts=find(diffindex>1);
                begin=indexexceed(pnts+1);
                ending=indexexceed(pnts);
                begin2=[indexexceed(1) begin];
                ending2=[ending indexexceed(end)];
                
                ind_b=[];
                for i=1:(length(begin2))
                    if (ending2(i)-begin2(i))>=100 % min duration of bursts
                        ind_b=[ind_b i];
                    end
                end
                
                if ~isempty (ind_b)
                    begin3=begin2(ind_b);
                    ending3=ending2(ind_b);
                    
                    space_betb=200; % min space between bursts
                    ind_b1=[];
                    for i=1:(length(begin3)-2)
                        if (begin3(i+1)-ending3(i))>=space_betb && (begin3(i+2)-ending3(i+1))>=space_betb
                            ind_b1=[ind_b1 i+1];
                        end
                    end
                    if (begin3(2)-ending(1))>=space_betb
                        ind_b1= [1 ind_b1];
                    end
                    if (begin3(length(begin3))-ending3(length(begin3)-1))>=space_betb
                        ind_b1=[ind_b1 length(begin3)];
                    end
                    onset1=begin3(ind_b1);
                    offset1=ending3(ind_b1);
                    %-------
                    epoch=200;
                end
                
                for j=1:length(onset1)
                    if onset1(j)>epoch && onset1(j)+epoch<length(env)
                        pre_b(j,1)=abs(mean(exp(sqrt(-1).*(dif_angs(rr,onset1(j):onset1(j)+epoch))))); %psi per burst
                        pre_ctx_b(j,:)=env(onset1(j)-epoch:onset1(j)+epoch);
                    end
                end
                
                for n=1:(length(env)/2)
                    idx_sur=randi([epoch+1,(length(env)-epoch)],1,1);
                    pre_sur(n,1)= abs(mean(exp(sqrt(-1).*(dif_angs(rr,idx_sur:idx_sur+epoch)))));
                end
                clust_b1(pp,1)=mean(pre_b); %psi per animal
                clust_s1(pp,1)=mean(pre_sur);
                ctx_b1(r,:)=mean(pre_ctx_b);
                pp=pp+1;
                
                all_psi_b=[all_psi_b; pre_b];
                all_psi_sur=[all_psi_sur; pre_sur];
                
                clust_b(r,:)=mean(clust_b1);
                clust_s(r,:)=mean(clust_s1);
                
                
            end
            
        end
        clearvars clust_b1 clust_s1 
        psi_fb{t,1}=all_psi_b;
        psi_fs{t,1}=all_psi_sur;
        
        if ~isempty (clust_b)
            
            psi_b(t,:)=mean(nonzeros(clust_b));%psi per frq
            psi_bsem(t,:)=std(nonzeros(clust_b))./sqrt(size((nonzeros(clust_b)),1));
            psi_s(t,:)=mean(nonzeros(clust_s));
            psi_ssem(t,:)=std(nonzeros(clust_s))./sqrt(size((nonzeros(clust_s)),1));
            ctx_b(t,:)=mean(nonzeros(ctx_b1));
            ctx_sb(t,:)=std(nonzeros(ctx_b1))./sqrt(size((nonzeros(ctx_b1)),1));
            
            stats(t,:)=ranksum(psi_fb{t,1},psi_fs{t,1})<= 0.05/size(psi_fb{t,1},1)
            
        else
            
            psi_b(t,:)=NaN;
            psi_s(t,:)=NaN;
            
        end
    end
    clearvars dif_angs pre_b pre_sur pre_ctx_b1
    ctcts_all(t,1)=ctcts;
    
end
freq_lab=find(ctcts_all~=0);
labels=({'5-15','16-26','27-37','38-48','49-60'})';
bar([psi_b psi_s])
legend('in bursts','surrogates')
xticklabels ({'5-15','16-26','27-37','38-48','49-60'});
% (labels{freq_lab,1})
xlabel ('Frequancy(Hz)')
ylabel ('PSI CTX-SNR during subctx bursts')
box off


% cd('/Users/Carolina/Desktop/TC_data') 
% savefig('subctxb_ctx_snr_01')
% saveas(gca,'subctxb_ctx_snr_01.png')

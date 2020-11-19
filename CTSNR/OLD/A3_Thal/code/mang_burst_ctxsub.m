clear all

cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
load ('data_all' , 'freq')
load 'SNR_ctx_probe'

samprate=1000;

for t=1:size(freq,1)
    clust_bi=[];
    clust_bo=[];
    all_ang_bi=[];
    all_ang_bo=[];
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
            %             plot(f,power)
            %             xlim ([0 100])
            
            
            [b,a]=butter(2,[(freq(t)-5)/(0.5*samprate) (freq(t)+5)/(0.5*samprate)],'bandpass');
            if t==length(freq)
                [b,a]=butter(2,[49/(0.5*samprate) 60/(0.5*samprate)],'bandpass');
            end
            [Pxx_ind,F_ind]=mscohere(data{r,1}(1,:),data{r,1}(rr,:),samprate,[],samprate,samprate); %Magnitude-squared coherence between ctx-a given contact
            m=1;
            if   sum(Pxx_ind(freq(t)-5:freq(t)+5))./sum(Pxx_ind(1:end))>0.1; %sum(power(freq(t)-5:freq(t)+5))./sum(power(1:end))>0.1
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
                        dumo=dif_angs(rr,onset1(j)-epoch:onset1(j));
                        dumi=dif_angs(rr,onset1(j):onset1(j)+(round((offset1(j)-onset1(j))/2)));
                        pre_bi(j,1)= (sum(exp(sqrt(-1)*(dumi))))./(size(dumi,2));
                        pre_bo(j,1)= (sum(exp(sqrt(-1)*(dumo))))./(size(dumo,2));
                    end
                end
                
                
                clust_bi1(pp,1)=mean(pre_bi); %psi per animal
                clust_bo1(pp,1)=mean(pre_bo);
                pp=pp+1;
                
                
                all_ang_bi=[all_ang_bi; pre_bi];
                all_ang_bo=[all_ang_bo; pre_bo];
                
                
                clust_bi(r,:)=mean(clust_bi1); %psi per animal
                clust_bo(r,:)=mean(clust_bo1);
                
            end
            
        end
        clearvars clust_b11 clust_bo1
    end
    
    ang_bi{t,1}=all_ang_bi;
    ang_bo{t,1}=all_ang_bo;
    
    if ~isempty (clust_bi)
        all_bi(t,:)=mean(nonzeros(clust_bi));%psi per frq
        all_bo(t,:)=mean(nonzeros(clust_bo));%psi per frq
        
    else
        
        all_bi(t,:)=NaN;
        all_bo(t,:)=NaN;
        
    end
    
    clearvars dif_angs pre_bi pre_bo
    ctcts_all(t,1)=ctcts;
    
end

angle(all_bi-all_bo)
for i=1:2
    for f=1:size(ang_bi,1)
        if i==1
            figure(i)
            subplot(size(ang_bi,1),1,f)
            polarhistogram(angle(ang_bi{f,1}))
        else
            figure(i)
            subplot(size(ang_bi,1),1,f)
            polarhistogram(angle(ang_bo{f,1}))
        end
    end
end



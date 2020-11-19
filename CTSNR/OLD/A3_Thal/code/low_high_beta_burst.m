clear all
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
% cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
A=char('BZ_ctx_probe.mat ','SNR_ctx_probe.mat');
B=char('bz_allrats.mat','snr_allrats.mat')
 rat_id=([1 7 9 11; 1 3 11 12]);
% rat_id=[([randi([1,5],1,1) randi([6,7],1,1) randi([8,10],1,1) 11]);([randi([1,2],1,1) randi([3,5],1,1) randi([6,11],1,1) randi([12,14],1,1)])];
samprate=1000;
colors = { [0.5 0 0.5] ; [0 0 0.5]};

% freq={[15:20];[21:30]};
freq={[15:30]};
for xx=1:size(A,1)
    cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
    name=A(xx,1:(find(A(xx,:)=='.')-1));
    load(name)
    cd ('C:\Users\creis\Documents\GitHub\CRcode\codes_thal\A4_Thal')
    name1=B(xx,1:(find(B(xx,:)=='.')-1));
    load(name1,'animals')
    color_b=colors{xx,1};
    data=data(animals(rat_id(xx,:)));
    clearvars -except xx rat_id colors A B freq samprate data color_b psi_all ctcts_all ctx_b psi_s psi_b stats
    for t=1:size(freq,1)
        ctx_b1=[];
        clust_b=[];
        clust_s=[];
        ctcts=0;
        for r=1:size(data,1) %electrodes
            thal_contact=[];
            coh_thal=[];
            filt_thal=[];
            clust_b1=[];
            clust_s1=[];
            clearvars a b
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
                    
                    dif_angs(rr,:)=non_norm;
                    clearvars non_norm;
                    
                    %                     channel_bursts=data{r,1}(rr,:);
                    channel_bursts=data{r,1}(1,:);
                    
                    [power1,f]=pwelch(channel_bursts,1000,[],1000,1000); %ctx power spectra
                    filt_bursts= find(power1==(max(power1(15:30,1))));
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
                    
                    if ranksum(pre_b,pre_sur)<0.05
                        clust_b1(pp,1)=mean(pre_b); %psi per animal
                        clust_s1(pp,1)=mean(pre_sur);
                        pp=pp+1;
                    end
                    ctx_b1(r,:)=mean(pre_ctx_b);
                    clust_b(r,:)=mean(clust_b1);
                    clust_s(r,:)=mean(clust_s1);
                end
            end
            if ~isempty (clust_b)
                psi_b(xx,t,:)=nanmean(nonzeros(clust_b));%psi per frq
                psi_s(xx,t,:)=nanmean(nonzeros(clust_s));
                ctx_b(xx,t,:)=mean(nonzeros(ctx_b1));
                stats(xx,t)=ranksum(nonzeros(clust_b(~isnan (clust_b))),nonzeros(clust_s(~isnan (clust_s))));
            else
                psi_b(xx,t,:)=NaN;
                psi_s(xx,t,:)=NaN;
                stats(xx,t)=NaN;
            end
        end
        clearvars dif_angs pre_b pre_sur pre_ctx_b1
        ctcts_all(xx,t)=ctcts;
    end
end

psi=[psi_b(1,:) psi_b(2,:);psi_s(1,:) psi_s(2,:)];
psi=[psi_b(1) psi_s(1);psi_b(2) psi_s(2)];

bar(psi)
xticklabels ({'BZ-CTX in {\beta}','SNr-CTX in {\beta}'})
ylabel ('Phase Synchrony Index')
box off
ylim ([0 1])

xticklabels ({'BZ-CTX','SNr-CTX',});
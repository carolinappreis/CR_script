clear all
close all
% cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
load ('data_all')
load 'bz_ctx_probe'
freq=[15;26];
color_bz=[0.5 0 0.5];
color_snr=[0 0 0.5];
color_ctxb=[0.5 0.5 0.5];

samprate=1000;

for t=1;
    %:size(freq,1)
    [b,a]=butter(2,[(freq(t)-5)/(0.5*samprate) (freq(t)+5)/(0.5*samprate)],'bandpass');
    ctx_b1=[];
    clust_b=[];
    clust_s=[];
    ctcts=0;
    for r=1:size(data,1); %electrodes
        thal_contact=[];
        coh_thal=[];
        filt_thal=[];
        m=1;
        [power1,f1]=pwelch(data{r,1}(1,:),1000,[],1000,1000);
        zs_power1=zscore(power1);
        for rr=2:size(data{r,1},1) %contacts
            [power,f]=pwelch(data{r,1}(rr,:),1000,[],1000,1000); %ctx power spectracontact power spectra
            zs_power=zscore(power);
            
            %             plot(f(1:50),log(power(1:50)))
            [Pxx_ind,F_ind]=mscohere(data{r,1}(1,:),data{r,1}(rr,:),samprate,[],samprate,samprate); %Magnitude-squared coherence between ctx-a given contact
            %             [(sum(Pxx_ind(15:30))./sum(Pxx_ind(1:end)))  sum(zs_power(15:30)>1.96)  sum(zs_power1(15:30)>1.96)]
            %             (sum(Pxx_ind(15:30))./sum(Pxx_ind(1:end))>0.1) && (sum(zs_power(15:30)>1.96)==0) && (sum(zs_power1(15:30)>1.96)==0)
            %             if  (sum(Pxx_ind(15:30))./sum(Pxx_ind(1:end))>0.1) && (sum(zs_power(15:30)>1.96)==0) && (sum(zs_power1(15:30)>1.96)==0); %sum(power(freq(t)-5:freq(t)+5))./sum(power(1:end))>0.1
            
            if  sum(Pxx_ind(15:30))./sum(Pxx_ind(1:end))>0.1
                thal_contact=[thal_contact rr];
                coh_thal=[coh_thal Pxx_ind(15:30)./sum(Pxx_ind(1:end))];
                filt_thal=[filt_thal ; filtfilt(b,a,data{r,1}(rr,:))]; %filt coherent subcortical contacts
                power_coh{r,1}(m,:)=power;
                power_ctx{r,1}(m,:)=power1;
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
                
                %                 channel_bursts=data{r,1}(rr,:);
                channel_bursts=data{r,1}(1,:);
                
                [power1,f]=pwelch(channel_bursts,1000,[],1000,1000); %ctx power spectra
                filt_bursts= find(power1==(max(power1(13:30,1))));
                [bb,aa]=butter(2,[(filt_bursts-5)/(0.5*samprate) (filt_bursts+5)/(0.5*samprate)],'bandpass');
                env=abs(hilbert(filtfilt(bb,aa,channel_bursts)));
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
                        pre_b(j,:)=dif_angs(rr,onset1(j)-epoch:onset1(j)+epoch);
                        pre_ctx_b(j,:)=env(onset1(j)-epoch:onset1(j)+epoch);
                    end
                end
                
                for n=1:(length(env)/100)
                    idx_sur=randi([epoch+1,(length(env)-epoch)],1,1);
                    pre_sur(n,:)= dif_angs(rr,idx_sur-epoch:idx_sur+epoch);
                end
                clust_b{r,1}(m,:)=abs(mean(exp(sqrt(-1).*(pre_b)))); %contact psi during bursts
                clust_s{r,1}(m,:)=abs(mean(exp(sqrt(-1).*(pre_sur))));
                ctx_b1(r,:)=mean(pre_ctx_b);
                m=m+1;
            end
        end
    end
    
    clust_b=clust_b(~cellfun('isempty',clust_b));
    animals=find(~cellfun('isempty',clust_b));
    clust_s=clust_s(~cellfun('isempty',clust_s));
    ctx_b1=ctx_b1(find(ctx_b1(:,1)~=0),:);
    power_coh= power_coh(~cellfun('isempty',power_coh));
    power_ctx= power_ctx(~cellfun('isempty',power_ctx));
    
%     for i =1:size(power_coh,1)
%         figure(1)
%         subplot(size(power_coh,1),1,i)
%         plot(log(mean(power_coh{i,1},1)))
%         figure(2)
%         subplot(size(power_ctx,1),1,i)
%         plot(log(mean(power_ctx{i,1},1)))%             plot(mean(clust_s{i,1},1))
%     end
    
    for cb=1:size(clust_b,1)
        psi_b(cb,:)=mean(clust_b{cb,1},1); %electrode psi (mean across contacts)
        psi_bsem(cb,:)=std(clust_b{cb,1})./sqrt(size(clust_b{cb,1},1));
        psi_s(cb,:)=mean(clust_s{cb,1},1);
        psi_ssem(cb,:)=std(clust_s{cb,1})./sqrt(size(clust_s{cb,1},1));
        
        ctx_b=mean(ctx_b1);
        ctx_sb=std(ctx_b1)./sqrt(size(ctx_b1,1));
        
        time= [1:2*epoch+1];
%         color_freq={[0 0 0.5],[0 0.5 0],[0.5 0 0],[0.8 0.2 0.8],[0.8 0.8 0]};
        subplot(size(clust_b,1)+1,1,1)
        y7=ctx_b; y6=y7+ctx_sb; y8=y7-ctx_sb;
%         color_b=color_freq{t}; color_s=[0 0 0];
        plot(time, y7, 'DisplayName','BZ aligned to short bursts')
        set(plot(time,y7),'LineWidth',1.5,'Color',color_ctxb);
        patch([time fliplr(time)], [y6 fliplr(y7)],[color_ctxb],'FaceAlpha',[0.2],'EdgeColor','none')
        patch([time fliplr(time)], [y7 fliplr(y8)],[color_ctxb],'FaceAlpha',[0.2],'EdgeColor','none')
        title ('PSI BZ-CTX at 10-20Hz in beta ctx bursts')
        xlim ([0 400])
        xticks([0:100:400])
        xticklabels ({'-200','-100','0','100','200'})
        hold on
        
        subplot(size(clust_b,1)+1,1,cb+1)
        y2=psi_b(cb,:); y1=y2+psi_bsem(cb,:); y3=y2-psi_bsem(cb,:);
        y5=psi_s(cb,:); y4=y5+psi_ssem(cb,:); y6=y5-psi_ssem(cb,:);
        
%         color_b=color_freq{t}; color_s=[0 0 0];
        color_b=color_bz; color_s=[0 0 0];
        plot(time, y2, 'DisplayName','BZ aligned to short bursts')
        set(plot(time,y2),'LineWidth',1.5,'Color',color_b);
        patch([time fliplr(time)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
        patch([time fliplr(time)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
        hold on
        plot(time, y5, 'DisplayName','BZ aligned to short bursts')
        set(plot(time,y5),'LineWidth',1.5,'Color',color_s);
        patch([time fliplr(time)], [y4 fliplr(y5)],[color_s],'FaceAlpha',[0.2],'EdgeColor','none')
        patch([time fliplr(time)], [y5 fliplr(y6)],[color_s],'FaceAlpha',[0.2],'EdgeColor','none')
        xlim ([0 400])
        ylim ([0 0.6])
        xticks([0:100:400])
        xticklabels ({'-200','-100','0','100','200'})
        
        
%         cd('C:\Users\creis\Documents\GitHub\CRcode\codes_thal\A4_Thal\code')
                 cd ('/Users/Carolina/Documents/GitHub/CRcode/codes_thal/A4_Thal/code')
        if size(clust_b{cb,1},1)>1
            st=NaN(1,401);
            clear A; A=clust_b{cb,:}; %b1{f,1};
            clear B; B=clust_s{cb,:}; %s1{f,1}(1:size(A,1),:);
            hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
            beg=find(st(1,:)<0.01 & st2(1,:)~=0);
            if ~isempty(beg)
                patch([beg(1) beg(end) beg(end) beg(1)],[0.55 0.55 0.575 0.575],color_b,'EdgeColor','none')
            end
        end
    end
    clearvars dif_angs a b
end


%
% cd('/Users/Carolina/Desktop/TC_data')
% savefig('clut_ctxb_ctx_snr')
% saveas(gca,'clust_ctxb_ctx_snr.png')

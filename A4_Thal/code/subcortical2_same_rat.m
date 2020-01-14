
clear all
close all
Fs=1000;
%  cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
 cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
load BZ_ctx_probe.mat
WaveA=WaveData_DCall;
dataA=data;
idx_A=[];
for i =1:29
    if size(WaveA{i,1},1)~=1
        idx_A=[idx_A i];
    end
end

%  cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
 cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
load SNr_ctx_probe.mat
WaveB=WaveData_DCall;
dataB=data;
idx_B=[];
for i =1:29
    if size(WaveB{i,1},1)~=1
        idx_B=[idx_B i];
    end
end

% load ('data_all.mat','freq')
freq=[5:35];
subj= idx_B(ismember(idx_B,idx_A));
samprate=1000;
ctx_b1=[];
clust_b=[];
clust_s=[];
dif_angs=[];
ctcts=0;
[b,a]=butter(2,[5/(0.5*samprate) 30/(0.5*samprate)],'bandpass');
for r=1:length(subj);
    m=1;
    powerA=[];
    powerB=[];
    for ra=2:size(WaveA{subj(r),1},1);
        for rb=2:size(WaveB{subj(r),1},1);
            [powerA,f]=pwelch(WaveA{subj(r),1}(ra,:),1000,[],1000,1000); %ctx power spectracontact power spectra
            [powerB,f]=pwelch(WaveB{subj(r),1}(rb,:),1000,[],1000,1000); %ctx power spectracontact power spectra
            [Pxx_indA,F_ind]=mscohere(WaveA{subj(r),1}(ra,:),WaveA{subj(r),1}(1,:),samprate,[],samprate,samprate); %Magnitude-squared coherence between ctx-a given contact
            [Pxx_indB,F_ind]=mscohere(WaveB{subj(r),1}(rb,:),WaveB{subj(r),1}(1,:),samprate,[],samprate,samprate); %Magnitude-squared coherence between ctx-a given contact
            
            
            %                 [sum(Pxx_indA(freq(t)-5:freq(t)+5))./sum(Pxx_indA(1:end)) sum(Pxx_indB(freq(t)-5:freq(t)+5))./sum(Pxx_indB(1:end))]
            %                 if  sum(Pxx_indA(freq(t)-5:freq(t)+5))./sum(Pxx_indA(1:end))>0.1 && sum(Pxx_indB(freq(t)-5:freq(t)+5))./sum(Pxx_indB(1:end))>0.1;
            if sum(Pxx_indA(freq))./sum(Pxx_indA(1:end))>0.1 && sum(Pxx_indB(freq))./sum(Pxx_indB(1:end))>0.1
                ctcts=ctcts+1;
                power_A{r,1}(m,:)=powerA;
                power_B{r,1}(m,:)=powerB;
                non_norm=squeeze(angle(hilbert(filtfilt(b,a,WaveA{subj(r),1}(ra,:))))-angle(hilbert(filtfilt(b,a,WaveB{subj(r),1}(rb,:))))); %dif angles between ctx-subctx
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
                
                  %                 channel_bursts=WaveB{r,1}(rb,:);
                %                 channel_bursts=WaveA{r,1}(ra,:);
                channel_bursts=WaveA{r,1}(1,:);
                
                freq_b=[15:35];
                [power1,f]=pwelch(channel_bursts,1000,[],1000,1000); %ctx power spectra
                filt_bursts= find(power1==(max(power1(freq_b,1))));
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
                    epoch=500;
%                     contacts_rat{r,1}=thal_contact;
                end
                
                for j=1:length(onset1)
                    if onset1(j)>epoch && onset1(j)+epoch<length(env)
                        pre_b(j,:)=dif_angs(onset1(j)-epoch:onset1(j)+epoch);
                        pre_ctx_b(j,:)=env(onset1(j)-epoch:onset1(j)+epoch);
                    end
                end
                
                for n=1:(length(env)/100)
                    idx_sur=randi([epoch+1,(length(env)-epoch)],1,1);
                    pre_sur(n,:)= dif_angs(idx_sur-epoch:idx_sur+epoch);
                end
                clust_b{r,1}(m,:)=abs(mean(exp(sqrt(-1).*(pre_b)))); %contact psi during bursts
                clust_s{r,1}(m,:)=abs(mean(exp(sqrt(-1).*(pre_sur))));
                ctx_b1(r,:)=mean(pre_ctx_b);
                m=m+1;
                dif_angs=[];
            end
        end
    end
end

clust_b=clust_b(~cellfun('isempty',clust_b));
animals=find(~cellfun('isempty',clust_b));
clust_s=clust_s(~cellfun('isempty',clust_s));
ctx_b1=ctx_b1(find(ctx_b1(:,1)~=0),:);
power_A= power_A(~cellfun('isempty',power_A));
power_B= power_B(~cellfun('isempty',power_B));

for cb=1:size(clust_b,1)
    psi_b(cb,:)=mean(clust_b{cb,1},1); %electrode psi (mean across contacts)
    psi_bsem(cb,:)=std(clust_b{cb,1})./sqrt(size(clust_b{cb,1},1));
    psi_s(cb,:)=mean(clust_s{cb,1},1);
    psi_ssem(cb,:)=std(clust_s{cb,1})./sqrt(size(clust_s{cb,1},1));
    
    ctx_b=mean(ctx_b1);
    ctx_sb=std(ctx_b1)./sqrt(size(ctx_b1,1));
end


color_b=[0.2 0.5 0.1];id_rat=animals; rats=size(id_rat,1);
color_s=[0 0 0];
color_ctxb=[0.5 0.5 0.5];
time= [1:2*epoch+1];
 
%-----

for i =1:length(id_rat)
clust_b{id_rat(i),:};
clust_b_m(i,:)=mean(ans);
clust_b_sd(i,:)=std(ans)./sqrt(size(ans,1));
clust_s{id_rat(i),:};
clust_s_m(i,:)=mean(ans);
clust_s_sd(i,:)=std(ans)./sqrt(size(ans,1));
power_B{id_rat(i),:};
power_B_m(i,:)=mean(ans);
power_B_sd(i,:)=std(ans)./sqrt(size(ans,1));
power_A1(i,:)=power_A{id_rat(i),1}(1,:);
end

ctx_b=ctx_b1(id_rat,:);
ctx_b_m=mean(ctx_b,1);
ctx_b_sd=std(ctx_b)./sqrt(size(ctx_b,1));

power_ctx_m=mean(power_A1);
power_ctx_sd=std(power_A1)./sqrt(size(power_A1,1));

%-----

reg_m=mean(clust_b_m);
reg_sd=std(clust_b_m)./sqrt(size(clust_b_m,1));
reg_sm=mean(clust_s_m);
reg_ssd=std(clust_s_m)./sqrt(size(clust_s_m,1));

cd('/Users/Carolina/Documents/GitHub/CR_script/A4_Thal/code') 
% cd('C:\Users\creis\Documents\GitHub\CRcode\codes_thal\A4_Thal\code')
st=NaN(1,2*epoch+1);
clear A; A=clust_b_m; %b1{f,1};
clear B; B=clust_s_m; %s1{f,1}(1:size(A,1),:);
hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
beg=find(st(1,:)<0.01 & st2(1,:)~=0);
if ~isempty(beg)
    sig_rise_all=[beg(1) beg(end)];
    
end

y2=reg_m; y1=y2+reg_sd; y3=y2-reg_sd;
y5=reg_sm; y4=y5+reg_ssd; y6=y5-reg_ssd;
p1=plot(time, y2,'LineWidth',1.5,'Color',color_b)
patch([time fliplr(time)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time fliplr(time)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
hold on
p2=plot(time, y5,'LineWidth',1.5,'Color',color_s);
patch([time fliplr(time)], [y4 fliplr(y5)],[color_s],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time fliplr(time)], [y5 fliplr(y6)],[color_s],'FaceAlpha',[0.2],'EdgeColor','none')
xlim ([0 400])
ylim ([0 0.6])
xticks([0:100:400])
xticklabels ({'-200','-100','0','100','200'})

patch([sig_rise_all(1) sig_rise_all(2) sig_rise_all(2) sig_rise_all(1)],[0 0 1 1],color_b,'FaceAlpha',[0.1],'EdgeColor','none')
title('BZ-CZ {\alpha} coupling during {\beta} bursts')
ylabel ('Phase Synchrony Index')
xlabel ('Time (msec)')
ylim([0 1])
box('off')
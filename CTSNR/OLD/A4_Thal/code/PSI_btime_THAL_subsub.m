
clear all
Fs=1000;
% cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
load SNR_ctx_probe.mat
WaveA=WaveData_DCall;
idx_A=[];
for i =1:29
    if size(WaveA{i,1},1)~=1
        idx_A=[idx_A i];
    end
end

% cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
load BZ_ctx_probe.mat
WaveB=WaveData_DCall;
idx_B=[];
for i =1:29
    if size(WaveB{i,1},1)~=1
        idx_B=[idx_B i];
    end
end

load ('data_all.mat','freq')
m=1;
subj= idx_B(ismember(idx_B,idx_A));
samprate=1000;

for t=1:size(freq,1)
    for r=1:length(subj)
        for n=2: size(WaveA{subj(r),1},1);
            [power,f]=pwelch(WaveA{subj(r),1}(n,:),1000,[],1000,1000);
            thal_contA(n,1)=power(freq(t));
        end
        sub_A=find(thal_contA==(max(thal_contA)));
        
        for o=2: size(WaveB{subj(r),1},1);
            [power,f]=pwelch(WaveB{subj(r),1}(o,:),1000,[],1000,1000);
            thal_contB(o,1)=power(freq(t));
        end
        sub_B=find(thal_contB==(max(thal_contB)));
        [b,a]=butter(2,[(freq(t)-5)/(0.5*Fs) (freq(t)+5)/(0.5*Fs)],'bandpass');
        if t==length(freq)
            [b,a]=butter(2,[49/(0.5*Fs) 60/(0.5*Fs)],'bandpass');
        end
        filtCTX(1,:)=filtfilt(b,a,WaveData_DCall{subj(r),1}(1,:));
        filtA(1,:)=filtfilt(b,a,WaveA{subj(r),1}(sub_A,:));
        filtB(1,:)=filtfilt(b,a,WaveB{subj(r),1}(sub_B,:));
        
        
        non_norm=squeeze(angle(hilbert(filtB))-angle(hilbert(filtA)));
        for x =1:size(non_norm,2)
            if non_norm(1,x)>pi
                non_norm(1,x)=(non_norm(1,x))-(2.*pi);
            elseif non_norm(1,x)<-pi
                non_norm(1,x)=(non_norm(1,x))+(2.*pi);
            else
                non_norm(1,x)= non_norm(1,x);
            end
        end
        dif_angs(r,:)=non_norm;
        psi_region=mean(abs(mean(exp(sqrt(-1)*(dif_angs)),2)));
        
        env=abs(hilbert(filtCTX));
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
            if (ending2(i)-begin2(i))>=100
                ind_b=[ind_b i];
            end
        end
        
        if ~isempty (ind_b)
            begin3=begin2(ind_b);
            ending3=ending2(ind_b);
            
            space_betb=200;
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
                pre_b(j,:)=dif_angs(r,onset1(j)-epoch:onset1(j)+epoch);
                pre_ctx_b(j,:)=env(onset1(j)-epoch:onset1(j)+epoch);
            end
        end
        
        pre_b( ~any( pre_b,2), : ) = [];
        pre_ctx_b( ~any( pre_ctx_b,2), : ) = [];
        
        for n=1:(length(env)/100)
            idx_sur=randi([epoch+1,(length(env)-epoch)],1,1);
            pre_sur(n,:)= dif_angs(r,idx_sur-epoch:idx_sur+epoch);
        end
        clust_b(r,:)=abs(mean(exp(sqrt(-1).*(pre_b))));
        clust_s(r,:)=abs(mean(exp(sqrt(-1).*(pre_sur))));
        ctx_b1(r,:)=mean(pre_ctx_b);
        
        clearvars thal_contA thal_contB
    end
    
        psi_b(t,:)=mean(clust_b,1);
        psi_bsem(t,:)=std(clust_b)./sqrt(size(clust_b,1));
        psi_s(t,:)=mean(clust_s,1);
        psi_ssem(t,:)=std(clust_s)./sqrt(size(clust_s,1));
        ctx_b(t,:)=mean(ctx_b1);
        ctx_sb(t,:)=std(ctx_b1)./sqrt(size(ctx_b1,1));
        
        time= [1:2*epoch+1];
        color_freq={[0 0 0.5],[0 0.5 0],[0.5 0 0],[0 0.5 0.5],[0.5 0.5 0]};
        subplot(6,1,1)
        y7=ctx_b(t,:); y6=y7+ctx_sb(t,:); y8=y7-ctx_sb(t,:);
        color_b=color_freq{t}; color_s=[0 0 0];
        plot(time, y7, 'DisplayName','BZ aligned to short bursts')
        set(plot(time,y7),'LineWidth',1.5,'Color',color_b);
        patch([time fliplr(time)], [y6 fliplr(y7)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
        patch([time fliplr(time)], [y7 fliplr(y8)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
        title ('PSI BZ-CZ in CTX bursts vs.surrogates')
        xlim ([0 400])
        xticks([0:100:400])
        xticklabels ({'-200','-100','0','100','200'})
        hold on
        
        subplot(6,1,t+1)
        y2=psi_b(t,:); y1=y2+psi_bsem(t,:); y3=y2-psi_bsem(t,:);
        y5=psi_s(t,:); y4=y5+psi_ssem(t,:); y6=y5-psi_ssem(t,:);
        
        color_b=color_freq{t}; color_s=[0 0 0];
        plot(time, y2, 'DisplayName','BZ aligned to short bursts')
        set(plot(time,y2),'LineWidth',1.5,'Color',color_b);
        patch([time fliplr(time)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
        patch([time fliplr(time)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
        hold on
        plot(time, y5, 'DisplayName','BZ aligned to short bursts')
        set(plot(time,y5),'LineWidth',1.5,'Color',color_s);
        patch([time fliplr(time)], [y4 fliplr(y5)],[color_s],'FaceAlpha',[0.2],'EdgeColor','none')
        patch([time fliplr(time)], [y5 fliplr(y6)],[color_s],'FaceAlpha',[0.2],'EdgeColor','none')
        %     ylim ([-0.05 0.3])
        %yticks ([0:0.1:0.2])
        xlim ([0 400])
        xticks([0:100:400])
        xticklabels ({'-200','-100','0','100','200'})
        
        cd('C:\Users\creis\Documents\GitHub\CRcode\codes_thal\A4_Thal\code')
        %   cd ('/Users/Carolina/Documents/GitHub/CRcode/codes_thal/A4_Thal/code')
        st=NaN(1,401);
        clear A; A=clust_b; %b1{f,1};
        clear B; B=clust_s; %s1{f,1}(1:size(A,1),:);
        hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
        beg=find(st(1,:)<0.05 & st2(1,:)~=0);
        if ~isempty(beg)
            patch([beg(1) beg(end) beg(end) beg(1)],[0.35 0.35 0.375 0.375],color_b,'EdgeColor','none')
        end
    clearvars dif_angs clust_b clust_s ctx_b1 
end


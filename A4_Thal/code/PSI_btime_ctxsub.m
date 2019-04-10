clear all
% cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
load ('data_all' , 'freq')
load 'AV_ctx_probe'

Fs=1000;
for t=1:size(freq,1)
    for ii=1:size(data,1)
        for r=2:size(data{ii,:},1)
            [b,a]=butter(2,[(freq(t)-5)/(0.5*Fs) (freq(t)+5)/(0.5*Fs)],'bandpass');
            if t==length(freq)
                [b,a]=butter(2,[49/(0.5*Fs) 60/(0.5*Fs)],'bandpass');
            end
            non_norm=angle(hilbert(filtfilt(b,a,data{ii,:}(1,:))))-(angle(hilbert(filtfilt(b,a,data{ii,:}(r,:)))));
            for x =1:size(non_norm,2)
                if non_norm(1,x)>pi
                    non_norm(1,x)=(non_norm(1,x))-(2.*pi);
                elseif non_norm(1,x)<-pi
                    non_norm(1,x)=(non_norm(1,x))+(2.*pi);
                else
                    non_norm(1,x)= non_norm(1,x);
                end
            end
            dif_angs1(r,:)=non_norm;
            clearvars non_norm;
        end
        dif_angs{ii,1}=dif_angs1(2:end,:);
        
        %---- bursts
        
        env=abs(hilbert(filtfilt(b,a,data{ii,:}(1,:))));
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
            
            for s=1:size(dif_angs{ii,:},1)
                for j=1:length(onset1)
                    if onset1(j)>epoch && onset1(j)+epoch<length(env)
                        pre_b{t,ii}(s,j,:)=dif_angs{ii,1}(s,onset1(j)-epoch:onset1(j)+epoch);
                        ctx_b1{t,ii}(s,j,:)=env(onset1(j)-epoch:onset1(j)+epoch);
                    end
                end
                
                for n=1:(length(env)/100)
                    idx_sur=randi([epoch+1,(length(env)-epoch)],1,1);
                    pre_sur{t,ii,1}(s,n,:)= dif_angs{ii,1}(s,idx_sur-epoch:idx_sur+epoch);
                end
            end
        end
        clust_b1{t,ii}=abs(mean(exp(sqrt(-1).*(reshape(pre_b{t,ii},[size(pre_b{t,ii},1)*size(pre_b{t,ii},2) size(pre_b{t,ii},3)])))));
        clust_s1{t,ii}=abs(mean(exp(sqrt(-1).*(reshape(pre_sur{t,ii},[size(pre_sur{t,ii},1)*size(pre_sur{t,ii},2) size(pre_sur{t,ii},3)])))));
        ctx_b1{t,ii}=mean(squeeze(mean(ctx_b1{t,ii},1)),1);
    end
    clust_b{t,1}= vertcat(clust_b1{t,:});
    clust_s{t,1}= vertcat(clust_s1{t,:});
    ctx_b(t,:)=mean(vertcat(ctx_b1{t,:}),1);
    ctx_sb(t,:)=std(vertcat(ctx_b1{t,:}),1)./sqrt(size(ctx_b1,1));
end

time= [1:2*epoch+1];
for f=1:size(pre_b,1)
    psi_b(f,:)=mean(clust_b{f,:},1);
    psi_bsem(f,:)=std(clust_b{f,:})./sqrt(size(clust_b{f,:},1));
    psi_s(f,:)=mean(clust_s{f,:},1);
    psi_ssem(f,:)=std(clust_s{f,:})./sqrt(size(clust_s{f,:},1));
end



for f=1:5
    color_freq={[0 0 0.5],[0 0.5 0],[0.5 0 0],[0 0.5 0.5],[0.5 0.5 0]};
    subplot(6,1,1)
    y7=ctx_b(f,:); y6=y7+ctx_sb(f,:); y8=y7-ctx_sb(f,:);
    color_b=color_freq{f}; color_s=[0 0 0];
    plot(time, y7, 'DisplayName','BZ aligned to short bursts')
    set(plot(time,y7),'LineWidth',1.5,'Color',color_b);
    patch([time fliplr(time)], [y6 fliplr(y7)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
    patch([time fliplr(time)], [y7 fliplr(y8)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
    xlim ([0 400])
    xticks([0:100:400])
    xticklabels ({'-200','-100','0','100','200'})
    title ('PSI CTX-AV in bursts vs.surrogates')
    hold on
    
    subplot(6,1,f+1)
    y2=psi_b(f,:); y1=y2+psi_bsem(f,:); y3=y2-psi_bsem(f,:);
    y5=mean(psi_s(f,:),1); y4=y5+psi_ssem(f,:); y6=y5-psi_ssem(f,:);
    
    color_b=color_freq{f}; color_s=[0 0 0];
    plot(time, y2, 'DisplayName','BZ aligned to short bursts')
    set(plot(time,y2),'LineWidth',1.5,'Color',color_b);
    patch([time fliplr(time)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
    patch([time fliplr(time)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
    hold on
    plot(time, y5, 'DisplayName','BZ aligned to short bursts')
    set(plot(time,y5),'LineWidth',1.5,'Color',color_s);
    patch([time fliplr(time)], [y4 fliplr(y5)],[color_s],'FaceAlpha',[0.2],'EdgeColor','none')
    patch([time fliplr(time)], [y5 fliplr(y6)],[color_s],'FaceAlpha',[0.2],'EdgeColor','none')
    ylim ([-0.05 0.3])
    %yticks ([0:0.1:0.2])
    xlim ([0 400])
    xticks([0:100:400])
    xticklabels ({'-200','-100','0','100','200'})
    
    cd('C:\Users\creis\Documents\GitHub\CRcode\codes_thal\A4_Thal\code')
    %   cd ('/Users/Carolina/Documents/GitHub/CRcode/codes_thal/A4_Thal/code')
    st=NaN(1,401);
    clear A; A=clust_b{f,1}; %b1{f,1};
    clear B; B=clust_s{f,1}; %s1{f,1}(1:size(A,1),:);
    hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
    beg=find(st(1,:)<0.05 & st2(1,:)~=0);
    if ~isempty(beg)
        patch([beg(1) beg(end) beg(end) beg(1)],[0.25 0.25 0.275 0.275],color_b,'EdgeColor','none')
    end
end





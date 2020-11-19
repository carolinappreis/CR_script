clear all
% cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
load ('data_all' , 'freq')
load 'CZ_ctx_probe'

Fs=1000;
for t=1:size(freq,1)
    for ii=1:size(data,1)
        for r=1:size(data{ii,:},1)
            [b,a]=butter(2,[(freq(t)-5)/(0.5*Fs) (freq(t)+5)/(0.5*Fs)],'bandpass');
            if t==length(freq)
                [b,a]=butter(2,[49/(0.5*Fs) 60/(0.5*Fs)],'bandpass');
            end
            
            filtall{ii,:}(r,:)=filtfilt(b,a,data{ii,:}(r,:));
            non_norm=squeeze(angle(hilbert(filtall{ii,1}(1,:)))-angle(hilbert(filtall{ii,1}(r,:))))';
            for x =1:size(non_norm,2)
                if non_norm(1,x)>pi
                    non_norm(1,x)=non_norm(1,x)-(2.*pi);
                elseif non_norm(1,x)<-pi
                    non_norm(1,x)=non_norm(1,x)+(2.*pi);
                end
            end
            dif_angs{ii,1}(r,:)=non_norm;
            
        end
        dif_angs{ii,1}=dif_angs{ii,1}(2:end,:);
        
        %---- bursts
        
        env=abs(hilbert(filtall{ii,1}(1,:)));
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
                if length(onset1)<epoch
                    for j=2:length(onset1)
                        pre_b{t,ii,1}(s,j-1)=(sum(exp(sqrt(-1)*(dif_angs{ii,1}(s,onset1(j)-epoch:onset1(j)+epoch))))./(length(onset1(j)-epoch:onset1(j)+epoch)));
                        %                     pre_b{t,ii,2}(s,j-1)=(sum(exp(sqrt(-1)*(dif_angs{t,ii,1}(s,onset1(j)-epoch:onset1(j)))))./(length(onset1(j)-epoch:onset1(j))));
                        
                    end
                elseif length(onset1)>epoch
                    for j=1:length(onset1)
                        pre_b{t,ii,1}(s,j)=(sum(exp(sqrt(-1)*(dif_angs{ii,1}(s,onset1(j)-epoch:onset1(j)+epoch))))./(length(onset1(j)-epoch:onset1(j)+epoch)));
                        %                     pre_b{t,ii,2}(s,j)=(sum(exp(sqrt(-1)*(dif_angs{t,ii,1}(s,onset1(j)-epoch:onset1(j)))))./(length(onset1(j)-epoch:onset1(j))));
                    end
                end
                
                for n=1:(length(env)/4)
                    idx_sur=randi([epoch+1,(length(env)-epoch)],1,1);
                    pre_sur{t,ii,1}(s,n,:)= (sum(exp(sqrt(-1)*(dif_angs{ii,1}(s,idx_sur:idx_sur+epoch)))))./(length(idx_sur:idx_sur+epoch));
                end
            end
        end
    end
end

for i =1:size(pre_b,1)
    for ii =1:size(pre_b,2)
        for j =1:size(pre_b{i,ii},1)
            pre_bm1(i,ii)=(mean(pre_b{i,ii}(:,:),1));   
        end
%         pre_surm(i,ii)=abs(mean(mean(pre_sur{i,ii}(:,:),2),1));
    end
end

bar([mean(psi_bm,2) mean(psi_surm,2)])
title ('PSI')
legend('in bursts','surrogates')
xticklabels({'5-15Hz','16-26Hz','27-37Hz','38-48Hz','49-100Hz'})
% 
% figure()
% for i=1:size(dif_angs,1)
%     n=[];
%     m=[];
%     for s=1:size(dif_angs,2)
%         n=[n double(ang_bm{i,s})];
%         subplot(1,size(ang_bm,1),i)
%         polarhistogram(n,12)
%     end
% end
% title ('mean angle inside burst')
% 
% figure()
% for i=1:size(dif_angs,1)
%     n=[];
%     m=[];
%     for s=1:size(dif_angs,2)
%         n=[n double(ang_surm{i,s})];
%         subplot(1,size(ang_bm,1),i)
%         polarhistogram(n,12)
%     end
% end
% title ('mean angle outside burst')






% ispc_trials=squeeze(abs(mean(exp(1i*diff(ang_trials,[],1)),2)));
% plot(1:1:size(output,3),ispc_trials)
%
% ispc_time= squeeze(abs(mean(exp(1i*diff(ang_trials,[],1)),3)));
% plot(1:1:size(output,2),ispc_time)


% cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A2_Thal/code')
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A2_Thal\code')
run 'psi_rat_BUA.m'

psi_ib=cell(2,1);
psi_ob=cell(2,1);
psi_sur=cell(2,1);
for rat =1:length(subj);
    
    env=abs(hilbert(filtCTX(rat,1,:)));
    threshold=prctile(env,75);
    tt(size(env,1):size(env,2))=threshold;
    
    indexexceed=find(env>threshold);
    diffindex=diff(indexexceed);
    pnts=find(diffindex>1);
    begin=indexexceed(pnts+1);
    ending=indexexceed(pnts);
    begin2=[indexexceed(1) begin'];
    ending2=[ending' indexexceed(end)];
    
    epoch=150;
    ind_b=[];
    for i=1:(length(begin2))
        if (ending2(i)-begin2(i))>=epoch
            ind_b=[ind_b i];
        end
    end
    
    begin3=begin2(ind_b);
    ending3=ending2(ind_b);
    
    duration=ending3-begin3;
    median_b=median(duration);
    SD_b=std(duration);
    
    space_betb=150;
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
    
    for s=1:size(dif_angs,2)
        for t=1:size(dif_angs,3)
            for j=1:length(onset1)
                psi_ib{rat,:}(s,t,j)=abs(sum(exp(sqrt(-1)*(dif_angs(rat,s,t,(onset1(j):onset1(j)+epoch)))))./(length(onset1(j):onset1(j)+epoch)));
                psi_ob{rat,:}(s,t,j)=abs(sum(exp(sqrt(-1)*(dif_angs(rat,s,t,(onset1(j)-epoch:onset1(j))))))./(length(onset1(j):onset1(j)+epoch)));
                pre_ang_ib{rat,:}(s,t,j)=(sum(exp(sqrt(-1)*(dif_angs(rat,s,t,(onset1(j):onset1(j)+epoch)))))./(length(onset1(j):onset1(j)+epoch)));
               pre_ang_ob{rat,:}(s,t,j)=(sum(exp(sqrt(-1)*(dif_angs(rat,s,t,(onset1(j)-epoch:onset1(j))))))./(length(onset1(j):onset1(j)+epoch)));
            end
        end
    end
    
    for s=1:size(dif_angs,2)
        for t=1:size(dif_angs,3)
            for n=1:length(env)./2
                idx_sur=randi([epoch+1,(length(env)-epoch)],1,1);
                psi_sur{rat,:}(s,t,n)= abs(sum(exp(sqrt(-1)*(dif_angs(rat,s,t,(idx_sur:idx_sur+epoch)))))./(length(idx_sur:idx_sur+epoch)));
                pre_ang_sur{rat,:}(s,t,n)=(sum(exp(sqrt(-1)*(dif_angs(rat,s,t,(idx_sur:idx_sur+epoch)))))./(length(idx_sur:idx_sur+epoch)));
            end
        end
    end
end

mean(mean(mean(psi_ib{f,:},3),2))
mean(mean(mean(psi_sur{f,:},3),2))
% mean(mean(mean(ang_ib{f,:},3),2))
% mean(mean(mean(ang_sur{f,:},3),2))

f=1;
imagesc(mean(psi_ib{f,:},3))
colorbar
figure()
imagesc(mean(psi_sur{f,:},3))
colorbar
figure()
imagesc(((mean(psi_ib{f,:}(:,:,:),3)-mean(psi_sur{f,:}(:,:,:),3))./(mean(psi_sur{f,:}(:,:,:),3))).*100)
colorbar
% bar(mean(psi_ib{1,:},3))
% figure()
% bar(mean(psi_sur{1,:},3))
% figure()
% bar(((mean(psi_ib{1,:}(:,:,:),3)-mean(psi_sur{1,:}(:,:,:),3))./(mean(psi_sur{1,:}(:,:,:),3))).*100)



% figure(1)
% for i=1:2
%     subplot(2,1,i)
% imagesc(squeeze(psi_ib{i,:}))
% imagesc(squeeze(psi_ib{i,:}))
% end
% figure(2)
% for i=1:2
%     subplot(2,1,i)
% imagesc(squeeze(psi_sur{i,:}))
% imagesc(squeeze(psisurr{i,:}))
% end
% 



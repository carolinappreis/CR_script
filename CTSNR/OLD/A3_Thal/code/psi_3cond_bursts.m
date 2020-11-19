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
    
    epoch=50;
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
    
     
   
    ind_b1=cell(2,1);
    for i=1:(length(begin3)-1)
        if (ending3(i)-begin3(i))<=median_b && (begin3(i+1)-ending3(i))>=200
            ind_b1{1}(1,i)=i;
        end
        if (ending3(i)-begin3(i))>=median_b && (begin3(i+1)-ending3(i))>=200
            ind_b1{2}(1,i)=i;
        end
    end
    
    space_betb=150;
    ind_b2=cell(2,1)
    for i=1:(length(begin3)-2)
        if (begin3(i+1)-ending3(i))>=space_betb && (begin3(i+2)-ending3(i+1))>=space_betb
            ind_b2=[ind_b2 i+1];
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
    
    
    for i=1:size(ind_b1,1)
        onset1{i,1}=begin3(ind_b1{i}(ind_b1{i}(1,:)~=0));
        offset1{i,1}=ending3(ind_b1{i}(ind_b1{i}(1,:)~=0));
    end
    
    duration1=cell(2,1);
    for i=1:2;
        duration1{i,1}=offset1{i}-onset1{i};
    end
end

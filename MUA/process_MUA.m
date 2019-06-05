
% cd('C:\Users\creis\Documents\GitHub\CR_script\MUA')
clear all
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\MUA')
load ('newMUA_SNR')
srn=1000;
spikerate=[];

for ii=1:size(MUA,1)
    for k=1:size(MUA{ii,1},1)
        sk_ctc=[];
        for i =1:srn:(length(MUA{ii,1}(k,:))-srn);
            sk_ctc=[sk_ctc numel(find(MUA{ii,1}(k,i:i+srn)==1))];
            SNR.spikerate{ii,1}(1,k)=mean(sk_ctc);
        end
        
        if SNR.spikerate{ii,1}(1,k)>=13 && SNR.spikerate{ii,1}(1,k)<=35
            beta_idx1{ii,1}(1,k)=k;
            b=beta_idx1{ii,1};
            SNR.beta_idx{ii,1}=(nonzeros(b))';
            clear b
        elseif  SNR.spikerate{ii,1}(1,k)>35 && SNR.spikerate{ii,1}(1,k)<=80
            gamma_idx1{ii,1}(1,k)=k;
            g=gamma_idx1{ii,1};
            SNR.gamma_idx{ii,1}= (nonzeros(g))';
            clear g
        elseif  SNR.spikerate{ii,1}(1,k)>7 && SNR.spikerate{ii,1}(1,k)<=12
            alpha_idx1{ii,1}(1,k)=k;
            a= alpha_idx1{ii,1};
            SNR.alpha_idx{ii,1}= (nonzeros(a))';
            clear a
        end
    end
    
end

aa=[];
bb=[];
gg=[];
for i=1:size(SNR.beta_idx,1)
    if ~isempty (SNR.alpha_idx{i,1})
        aa=[aa;i];
    end
    
    if ~isempty (SNR.beta_idx{i,1})
        bb=[bb;i];
    end
    
    if ~isempty (SNR.gamma_idx{i,1})
        gg=[gg;i];
    end
end
SNR.ctx=ctx;
SNR.mua=MUA;
SNR.alpha_rats=aa;
SNR.beta_rats=bb;
SNR.gamma_rats=gg;

clearvars -except SNR
%cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\MUA\probe MUA_act_mat')
% save ('SNR_spikerate.mat')

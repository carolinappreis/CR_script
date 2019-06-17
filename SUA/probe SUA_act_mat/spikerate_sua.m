clear all
% cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/probe SUA_act_mat')
cd ('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\probe SUA_act_mat')
load('probe_SUA_CZ.mat')

srn=1000;
spikerate=[];

SUA=data_all; 
ctx=Ecog_all;

for ii=1:size(SUA,1)
    for k=1:size(SUA{ii,1},1)
        sk_ctc=[];
        for i =1:srn:(length(SUA{ii,1}(k,:))-srn);
            sk_ctc=[sk_ctc numel(find(SUA{ii,1}(k,i:i+srn)==1))];
            CZ.spikerate{ii,1}(1,k)=mean(sk_ctc);
        end
        
%         if CZ.spikerate{ii,1}(1,k)>=13 && CZ.spikerate{ii,1}(1,k)<=35
%             beta_idx1{ii,1}(1,k)=k;
%             b=beta_idx1{ii,1};
%             CZ.beta_idx{ii,1}=(nonzeros(b))';
%             clear b
%         elseif CZ.spikerate{ii,1}(1,k)>35 && CZ.spikerate{ii,1}(1,k)<=80
%             gamma_idx1{ii,1}(1,k)=k;
%             g=gamma_idx1{ii,1};
%             CZ.gamma_idx{ii,1}= (nonzeros(g))';
%             clear g
%         elseif  CZ.spikerate{ii,1}(1,k)>7 && CZ.spikerate{ii,1}(1,k)<=12
%             alpha_idx1{ii,1}(1,k)=k;
%             a= alpha_idx1{ii,1};
%             CZ.alpha_idx{ii,1}= (nonzeros(a))';
%             clear a
%         end
    end
    
end

aa=[];
bb=[];
gg=[];
for i=1:size(CZ.alpha_idx,1)
    if ~isempty (CZ.alpha_idx{i,1})
        aa=[aa;i];
    end
end

for i=1:size(CZ.beta_idx,1)
    if ~isempty (CZ.beta_idx{i,1})
        bb=[bb;i];
    end
end

% for i=1:size(CZ.gamma_idx,1)
%     if ~isempty (CZ.gamma_idx{i,1})
%         gg=[gg;i];
%     end
% end

CZ.ctx=ctx;
CZ.sua=SUA;
CZ.alpha_rats=aa;
CZ.beta_rats=bb;
% CZ.gamma_rats=gg;

clearvars -except CZ
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\probe SUA_act_mat')
save ('CZ_sua_skrate.mat')

% clear all
cd ('C:\Users\wolf5173\Documents\Carolina_code\codes_thal\A2_Thal\mat')
load 'animal_lesion_nolesion.mat'

beta_coherent=beta_coh_all;

for j=1:size(beta_coh_all,1)
    if ~isempty (beta_coh_all{j});
        for k=1:length(beta_coh_all{j,1})
            if beta_coherent{j}(1,k)<0.1;
                beta_coherent{j}(1,k)=NaN;
            end
        end
    end
end

for j=1:size(beta_coh_all,1)
    if ~isempty (beta_coh_all{j});
        index_beta{j,:}=find (~isnan(beta_coherent{j}(1,:)));
    end
end

f=[];
for j=1:size(evoked_all,1)
if ~isempty (evoked_all{j})
f = [f j];
end
end


%we know that 'kjx140'  f01@70-170_m.mat and e01@70-170_m.mat' are the
%exact same location - decide for only one to stay

data= A(lesion(f),:);


% if size(index_beta{f(6)},2) < size(index_beta{f(7)},2)
%     
%     power_all {f(6)} = [];
%     evoked_all {f(6)}=[];
%     power_all {f(6)}=[];
%     beta_coh_all {f(6)}=[];
%     coh_all {f(6)}=[];
% else
%     power_all {f(7)} = [];
%     evoked_all {f(7)}=[];
%     power_all {f(7)}=[];
%     beta_coh_all {f(7)}=[];
%     coh_all {f(7)}=[];
% end
% 
% clearvars -except power_all evoked_all power_all beta_coh_all coh_all 
% save 'lesion_param_2000_BZ_new'

counting=[];
m=1;
for i =1:size(index_beta,1)
    if ~isempty(index_beta{i})
        counting(1,m)=size(index_beta{i},2);
        m=m+1;
    end
end
counting2=sum(counting)
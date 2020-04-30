clear all
close all

axi(1,:) = [NaN 3 1 1 NaN 1 3 NaN 1 2]; %% as ploted in pca plots
axi(2,:) = [1 1 3 1 1 NaN 1 1 1 1];
% main=[1 1 3 1 3 3 3 3 1 1];
ns_mat=[[1 2 3]; [1 2 3]; [3 2 1]; [1 2 3];[3 2 1]; [3 2 1]; [3 2 1]; [3 2 1]; [1 2 3]; [1 2 3]];

load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/clusters_BA.mat')
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/output_auto.mat','TT1_C1')
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/output_auto.mat','TT1_C2')

for numb=1:length(cohort)
    r(numb,:)=[numel(find(C_RS{numb,1}(:)==1));numel(find(C_RS{numb,1}(:)==2))]';
    rs_max(numb,1)=find(r(numb,:)==(max(r(numb,:))));
    
    n(numb,:)=[numel(find(C_NS{numb,1}(:)==1));numel(find(C_NS{numb,1}(:)==2))]';
    ns_max(numb,1)=find(n(numb,:)==(max(n(numb,:))));
    
    ttd=eval(['cat(3,TT1_C' num2str(rs_max(numb,1)) '{' num2str(numb) ',:})']);
    tt_all{numb,1}=ttd;
    tt1(numb,:,:)=(squeeze(nanmedian(ttd,1)))'; clear ttd;
    dumm= squeeze(eval(['nostim_c' num2str(rs_max(numb,1)) '(' num2str(numb) ',:,:)']));
    %     nostim(numb,:,:)=dumm(squeeze(ns_mat(rs_max(numb),numb,:)),:); clear dumm;
    nostim(numb,:,:)=dumm(ns_mat(numb,:),:); clear dumm;
    main_clust(numb)=eval(['axi(' num2str(rs_max(numb,1)) ',' num2str(numb) ')']);
    dumi=squeeze(eval(['nostimout_c' num2str(rs_max(numb,1)) '(' num2str(numb) ',:,:)']));
    nostimout(numb,:,:)=dumi(ns_mat(numb,:),:); clear dumi;
end

clearvars -except tt1 nostim nostimout main_clust tt_all
%output name data_matchcluster.mat
%raw ARC with NS thereshold for 3 axis








% if c == 1
%     %     if strcmp(method, 'ward-pca2')
%     axi = [NaN 3 1 1 NaN 1 3 NaN 1 2]; %% as ploted in pca plots
%     main=[NaN 3 3 1 NaN 3 1 NaN 1 2]; %% matching to ns stim (as connected in the beginning)
%     main=main(y);
%     ns_mat = [repmat(NaN,1,3); [3 2 1]; [3 2 1]; [1 2 3];repmat(NaN,1,3); [3 2 1]; [1 2 3]; repmat(NaN,1,3); [1 2 3]; [2 3 1]];
%     ns_mat=ns_mat(y,:);
%     %     end
% end
%
% if c == 2 %% MUST ADD CASE TO REMOVE FAULTY BASELINE
%     %     if strcmp(method, 'ward-pca2')
%     axi = [1 1 3 1 1 NaN 1 1 1 1];
%     main = [1 1 1 1 3 NaN 3 3 1 1];
%     main=main(y);
%     ns_mat = [[1 2 3]; [1 2 3]; [1 2 3]; [1 2 3]; [3 2 1]; repmat(NaN,1,3); [3 2 1]; [3 2 1]; [1 2 3]; [1 2 3]];
%     ns_mat=ns_mat(y,:);
%     %     end
% end
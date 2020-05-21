clear all
close all

axi(1,:) = [NaN 3 1 1 NaN 1 3 NaN 1 2]; %% as ploted in pca plots
axi(2,:) = [1 1 3 1 1 NaN 1 1 1 1];
% main=[1 1 3 1 3 3 3 3 1 1];
ns_mat=[[1 2 3]; [1 2 3]; [3 2 1]; [1 2 3];[3 2 1]; [3 2 1]; [3 2 1]; [3 2 1]; [1 2 3]; [1 2 3]];

load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/clusters_CR.mat')
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/output_auto.mat','TT1_C1')
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/output_auto.mat','TT1_C2')
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/segments_output.mat','TT1_C2','TT1_C1','TT1_A2','TT1_A1','seg_filt1','seg_filt2','seg_zfilt1','seg_zfilt2')
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/segments_output2.mat','TT1_C2','TT1_C1')

cohort = [ 2 3 4 5 8 10 11 13 16 17];
for numb=1:length(cohort)
    r(numb,:)=[numel(find(C_RS{numb,1}(:)==1));numel(find(C_RS{numb,1}(:)==2))]';
    
    if r(numb,1)==r(numb,2)
        rs_max(numb,1)=1;
    else
        rs_max(numb,1)=find(r(numb,:)==(max(r(numb,:)))); 
    end
    
    n(numb,:)=[numel(find(C_NS{numb,1}(:)==1));numel(find(C_NS{numb,1}(:)==2))]';
    ns_max(numb,1)=find(n(numb,:)==(max(n(numb,:))));
    
    ttd=eval(['cat(3,TT1_C' num2str(rs_max(numb,1)) '{' num2str(numb) ',:})']);
    tt_all{numb,1}=ttd;
    tt1(numb,:,:)=(squeeze(nanmedian(ttd,1)))'; clear ttd;
    
    %     ttd=eval(['cat(3,TT1_A' num2str(rs_max(numb,1)) '{' num2str(numb) ',:})']);
    %     tt_all_amp{numb,1}=ttd;
    %     tt1_a(numb,:,:)=(squeeze(nanmedian(ttd,1)))'; clear ttd;
    
    dumm= squeeze(eval(['nostim_c' num2str(rs_max(numb,1)) '(' num2str(numb) ',:,:)']));
    %     nostim(numb,:,:)=dumm(squeeze(ns_mat(rs_max(numb),numb,:)),:); clear dumm;
    nostim(numb,:,:)=dumm(ns_mat(numb,:),:); clear dumm;
%     main_clust(numb)=eval(['axi(' num2str(rs_max(numb,1)) ',' num2str(numb) ')']);
    %     dumi=squeeze(eval(['nostimout_c' num2str(rs_max(numb,1)) '(' num2str(numb) ',:,:)']));
    %     nostimout(numb,:,:)=dumi(ns_mat(numb,:),:); clear dumi;
    
    %     for ii=1:12
    %         segev{numb,ii}=squeeze(eval(['seg_filt' num2str(rs_max(numb,1)) '{' num2str(numb) ',' num2str(ii) '}(' num2str(main_clust(numb)) ',:,:)']));
    %         segevz{numb,ii}=squeeze(eval(['seg_zfilt' num2str(rs_max(numb,1)) '{' num2str(numb) ',' num2str(ii) '}(' num2str(main_clust(numb)) ',:,:)']));
    %     end
end

clearvars -except tt1 nostim nostimout main_clust tt_all segev tt_all_amp tt1_a r rs_max ns_max n 
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
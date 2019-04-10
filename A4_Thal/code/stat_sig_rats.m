
clear all
close all
run 'psi_rats_savemat.m'
clearvars -except stat_idx clust_b clust_s clust_s power_coh power_ctx ctx_b1 psi_b  psi_bsem psi_s psi_ssem ctx_b ctx_sb time 

for i =1:length(stat_idx)
clust_b{stat_idx(i),:};
clust_b_m(i,:)=mean(ans);
clust_b_sd(i,:)=std(ans)./sqrt(size(ans,1));
clust_s{stat_idx(i),:};
clust_s_m(i,:)=mean(ans);
clust_s_sd(i,:)=std(ans)./sqrt(size(ans,1));
power_coh{stat_idx(i),:};
power_coh_m(i,:)=mean(ans);
power_coh_sd(i,:)=std(ans)./sqrt(size(ans,1));
power_ctx1(i,:)=power_ctx{stat_idx(i),1}(1,:);
end

ctx_b=ctx_b1(stat_idx,:);
ctx_b_m=mean(ctx_b,1);
ctx_b_sd=std(ctx_b)./sqrt(size(ctx_b,1));

power_ctx_m=mean(power_ctx1);
power_ctx_sd=std(power_ctx1)./sqrt(size(power_ctx1,1));

% cd ('/Users/Carolina/Documents/GitHub/CRcode/codes_thal/A4_Thal/code')
for r=1:size(stat_idx,2);
    st=NaN(1,401);
    clear A; A=clust_b{stat_idx(r),:}; %b1{f,1};
    clear B; B=clust_s{stat_idx(r),:}; %s1{f,1}(1:size(A,1),:);
    hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
    beg=find(st(1,:)<0.01 & st2(1,:)~=0);
    if ~isempty(beg)
        sig_rise(r,:)=[beg(1) beg(end)];
    end
end


clearvars -except sig_rise stat_idx clust_b clust_b_m clust_b_sd clust_s clust_s_m clust_s_sd power_coh power_coh_m power_coh_sd power_ctx1 ctx_b ctx_b_sd ctx_b_m power_ctx_m power_ctx_sd ctx_b_sd
% % cd('C:\Users\creis\Documents\GitHub\CRcode\codes_thal\A4_Thal')
% cd('/Users/Carolina/Documents/GitHub/CRcode/codes_thal/A4_Thal')
% save ('snr_rat_level.mat')



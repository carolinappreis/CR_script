clear all; close all
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/aic_r_results.mat')
% aicx={'aic_c';'aic_n'};
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/aic_r4.mat')
% aicx={'Aic_raw4';'Aic_n4'};
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/aic_normalised_quartic.mat')
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/r_poly.mat')
aicx={'Aic_n4'};


for i=1
data=eval(aicx{i,1});
for m=1:size(data,1)
  [s_val s_idx]=sort(data(m,:),'ascend');
  win(m,i,1)=s_idx(1);
%   dif_aic(m,:)=[(abs(s_val(1)-s_val(2)))>2];
    dif_aic(m,:)=[(abs(s_val(1)-s_val(2)))>2 (abs(s_val(1)-s_val(3)))>2 (abs(s_val(1)-s_val(4)))>2];
 
end

end





% % 
% % clear all
% % 
% % load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cleaned_rc12_noaddon.mat')
% % load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/newnonstim10.mat')
% % 
% % for i=1:10
% %     s=nanmedian(tt1{i,1});
% %     ns=squeeze(nostimout(i,1,:))';
% %     st=NaN(1,12);
% %     A=s;
% %     B=ns;
% %     hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
% %     beg=find(st(1,:)<0.05 & st2(1,:)~=0);
% %     if ~isempty(beg)
% %         sig_rise_all=[beg(1) beg(end)];
% %         
% %     end
% %     
% %     
% % end
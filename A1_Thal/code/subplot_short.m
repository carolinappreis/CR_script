% Ploting median change in beta activity at the subcortical level alignedwith short bursts (beta coherent and non-coherent probes) 

clear all
% cd ('\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A1_Thal\mat')
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A1_Thal/mat')
figure1 = figure('InvertHardcopy','off','Color',[1 1 1],...
'OuterPosition',[1009 27.6666666666667 373.333333333333 843.333333333333]);
load ('bursts_lesioned.mat','bursts_plot')
ax8=subplot(8,1,1)
x=[-500:500];
y1=bursts_plot{1}(1,:);
y2=bursts_plot{1}(2,:);
y3=bursts_plot{1}(3,:);
color_BZ=[0.1 0.1 0.1];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_BZ);
patch([x fliplr(x)], [y1 fliplr(y2)],color_BZ,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_BZ,'FaceAlpha',0.2,'EdgeColor','none')


load 'BZ_nc_fig.mat'
ax1=subplot(8,1,2)
x=[-500:500];
y1=med_thal_s{1}(1,:);
y2=med_thal_s{1}(2,:);
y3=med_thal_s{1}(3,:);
color_BZ=[0.5 0 0.5];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_BZ);
patch([x fliplr(x)], [y1 fliplr(y2)],color_BZ,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_BZ,'FaceAlpha',0.2,'EdgeColor','none')
hold on
out3_real= all_contacts;
clear all_contacts
out3_sur= surr2;
clear surr2
color_BZ=[0.5 0 0.5];
color_SNR=[0 0 0.5];
color_region=color_BZ;
st=NaN(2,551);
for ind=1:2; clear A; A(1:size(out3_real,2),1:551)=out3_real(ind,:,300:850);
    clear B; B(1:size(out3_sur,1),1:551)=out3_sur(:,300:850);
    hayriye_c; st(ind,:)=stats.prob; st2(ind,:)=stats.posclusterslabelmat;
end

beg=find(st(1,:)<0.05 & st2(1,:)~=0);
if ~isempty(beg)
    beg(1)
    beg(find(diff(beg)>1))
    beg(find(diff(beg)>1)+1)
    beg(end)
    patch([beg(1) beg(end) beg(end) beg(1)]-200,[0.15 0.15 0.175 0.175],color_region,'EdgeColor','none')
end
hold on
x=[-500:500];
y1=surr_s(1,:);
y2=surr_s(2,:);
y3=surr_s(3,:);
color_BZ=[0 0 0];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_BZ);
patch([x fliplr(x)], [y1 fliplr(y2)],color_BZ,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_BZ,'FaceAlpha',0.2,'EdgeColor','none')
title ('BZ')


load 'SNr_nc_fig.mat'
ax2=subplot(8,1,3)
x=[-500:500];
y1=med_thal_s{1}(1,:);
y2=med_thal_s{1}(2,:);
y3=med_thal_s{1}(3,:);
color_SNR=[0 0 0.5];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_SNR);
patch([x fliplr(x)], [y1 fliplr(y2)],color_SNR,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_SNR,'FaceAlpha',0.2,'EdgeColor','none')
hold on

load ('SNr_nc_fig','all_contacts','surr2','surr_s')
out3_real= all_contacts;
clear all_contacts
out3_sur= surr2;
clear surr2
color_BZ=[0.5 0 0.5];
color_SNR=[0 0 0.5];
color_region=color_SNR;

st=NaN(2,551);
for ind=1:2; clear A; A(1:size(out3_real,2),1:551)=out3_real(ind,:,300:850);
    clear B; B(1:size(out3_sur,1),1:551)=out3_sur(:,300:850);
    hayriye_c; st(ind,:)=stats.prob; st2(ind,:)=stats.posclusterslabelmat;
end

beg=find(st(1,:)<0.05 & st2(1,:)~=0);
if ~isempty(beg)
    beg(1)
    beg(find(diff(beg)>1))
    beg(find(diff(beg)>1)+1)
    beg(end)
    patch([beg(1) beg(end) beg(end) beg(1)]-200,[0.15 0.15 0.175 0.175],color_region,'EdgeColor','none')
end
load 'SNr_nc_fig.mat'
x=[-500:500];
y1=surr_s(1,:);
y2=surr_s(2,:);
y3=surr_s(3,:);
color_SNR=[0 0 0];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_SNR);
patch([x fliplr(x)], [y1 fliplr(y2)],color_SNR,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_SNR,'FaceAlpha',0.2,'EdgeColor','none')
title ('SNr')



load 'CZ_nc_fig.mat'
ax3=subplot(8,1,4)
x=[-500:500];
y1=med_thal_s{1}(1,:);
y2=med_thal_s{1}(2,:);
y3=med_thal_s{1}(3,:);
color_CZ=[0 0.2 0.5];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_CZ);
patch([x fliplr(x)], [y1 fliplr(y2)],color_CZ,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_CZ,'FaceAlpha',0.2,'EdgeColor','none')
hold on
load 'CZ_nc_fig.mat'
x=[-500:500];
y1=surr_s(1,:);
y2=surr_s(2,:);
y3=surr_s(3,:);
color_CZ=[0 0 0];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_CZ);
patch([x fliplr(x)], [y1 fliplr(y2)],color_CZ,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_CZ,'FaceAlpha',0.2,'EdgeColor','none')
title ('CZ')

load 'AV_nc_fig.mat'
ax4=subplot(8,1,5)
x=[-500:500];
y1=med_thal_s{1}(1,:);
y2=med_thal_s{1}(2,:);
y3=med_thal_s{1}(3,:);
color_AV=[0 0.5 0];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_AV);
patch([x fliplr(x)], [y1 fliplr(y2)],color_AV,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_AV,'FaceAlpha',0.2,'EdgeColor','none')
hold on
load 'AV_nc_fig.mat'
x=[-500:500];
y1=surr_s(1,:);
y2=surr_s(2,:);
y3=surr_s(3,:);
color_AV=[0 0 0];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_AV);
patch([x fliplr(x)], [y1 fliplr(y2)],color_AV,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_AV,'FaceAlpha',0.2,'EdgeColor','none')
title ('AV')

% clear all
% load 'nc_AM_200p.mat'
% x=[-500:500];
% y1=med_thal_s{1}(1,:);
% y2=med_thal_s{1}(2,:);
% y3=med_thal_s{1}(3,:);
% color_AM=[0.5 0 0.5];
% plot(x,y2)
% set(plot(x,y2),'LineWidth',1.5,'Color',color_AM);
% patch([x fliplr(x)], [y1 fliplr(y2)],color_AM,'FaceAlpha',0.2,'EdgeColor','none')
% patch([x fliplr(x)], [y2 fliplr(y3)],color_AM,'FaceAlpha',0.2,'EdgeColor','none')
% 

load 'RT_nc_fig.mat'
ax5=subplot(8,1,6)
x=[-500:500];
y1=med_thal_s{1}(1,:);
y2=med_thal_s{1}(2,:);
y3=med_thal_s{1}(3,:);
color_RT=[0.2 0.7 0.5];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_RT);
patch([x fliplr(x)], [y1 fliplr(y2)],color_RT,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_RT,'FaceAlpha',0.2,'EdgeColor','none')
hold on
load 'RT_nc_fig.mat'
x=[-500:500];
y1=surr_s(1,:);
y2=surr_s(2,:);
y3=surr_s(3,:);
color_RT=[0 0 0];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_RT);
patch([x fliplr(x)], [y1 fliplr(y2)],color_RT,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_RT,'FaceAlpha',0.2,'EdgeColor','none')
title ('Rt')


% clear all
% load 'NW_nc_fig.mat'
% subplot(8,1,6)
% x=[-500:500];
% y1=med_thal_s{1}(1,:);
% y2=med_thal_s{1}(2,:);
% y3=med_thal_s{1}(3,:);
% color_NW=[0 0.1 0.2];
% plot(x,y2)
% set(plot(x,y2),'LineWidth',1.5,'Color',color_NW);
% patch([x fliplr(x)], [y1 fliplr(y2)],color_NW,'FaceAlpha',0.2,'EdgeColor','none')
% patch([x fliplr(x)], [y2 fliplr(y3)],color_NW,'FaceAlpha',0.2,'EdgeColor','none')
% title ('NW')

load 'PC_nc_fig.mat'
ax6=subplot(8,1,7)
x=[-500:500];
y1=med_thal_s{1}(1,:);
y2=med_thal_s{1}(2,:);
y3=med_thal_s{1}(3,:);
color_PC=[0.4 0 0];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_PC);
patch([x fliplr(x)], [y1 fliplr(y2)],color_PC,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_PC,'FaceAlpha',0.2,'EdgeColor','none')
hold on
load 'PC_nc_fig.mat'
x=[-500:500];
y1=surr_s(1,:);
y2=surr_s(2,:);
y3=surr_s(3,:);
color_PC=[0 0 0];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_PC);
patch([x fliplr(x)], [y1 fliplr(y2)],color_PC,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_PC,'FaceAlpha',0.2,'EdgeColor','none')
title ('PC')

load 'ZI_nc_fig.mat'
ax7=subplot(8,1,8)
x=[-500:500];
y1=med_thal_s{1}(1,:);
y2=med_thal_s{1}(2,:);
y3=med_thal_s{1}(3,:);
color_ZI=[0.75 0 0];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_ZI);
patch([x fliplr(x)], [y1 fliplr(y2)],color_ZI,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_ZI,'FaceAlpha',0.2,'EdgeColor','none')
xlabel('msec','FontSize',12)
hold on
load 'ZI_nc_fig.mat'
x=[-500:500];
y1=surr_s(1,:);
y2=surr_s(2,:);
y3=surr_s(3,:);
color_ZI=[0 0 0];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_ZI);
patch([x fliplr(x)], [y1 fliplr(y2)],color_ZI,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_ZI,'FaceAlpha',0.2,'EdgeColor','none')
title ('ZI')
xlabel('msec','FontSize',14)
ylabel('Change in beta power (%Baseline)','FontSize',14)

ylim([ ax1 ax2 ax3 ax4 ax5 ax6 ax7],[-0.2 0.2])
ylim(ax8, [-2 2])
xticks([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8],[-500:250:500])
yticklabels([ax1 ax2 ax3 ax4 ax5 ax6 ax7],[-20 0 20])



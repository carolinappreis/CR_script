 clear all
% cd ('/Users/Carolina/Documents/MATLAB/hayriye code for tc analysis')
% load('stn_amp_real','out3');
% out3_real=out3;
% clear out3
% load('stn_amp_sur_comb','out3');
% out3_sur=out3;
% clear out3

cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A1_Thal/mat')
load ('BZ_fig','med_contacts','surr','surr_s')
out3_real= med_contacts;
clear med_contacts
out3_sur= surr;
clear surr
color_BZ=[0.5 0 0.5];
color_SNR=[0 0 0.5];
color_surr=[0.3 0.3 0.3];
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
    figure(1)
    ax1=subplot(2,2,1)
    patch([beg(1) beg(end) beg(end) beg(1)]-200,[0.2 0.2 0.225 0.225],color_region,'EdgeColor','none')
end
hold on
clear beg
beg=find(st(2,:)<0.05 & st2(2,:)~=0);
if ~isempty(beg)
    beg(1)
    beg(find(diff(beg)>1))
    beg(find(diff(beg)>1)+1)
    beg(end)
    ax2=subplot(2,2,3)
    patch([beg(1) beg(end) beg(end) beg(1)]-200,[0.2 0.2 0.225 0.225],color_region,'EdgeColor','none')
end
hold on
load ('BZ_fig','med_thal_s')
figure(1)
ax1=subplot(2,2,1)
x=[-500:500];
y1=med_thal_s{1}(1,:);
y2=med_thal_s{1}(2,:);
y3=med_thal_s{1}(3,:);
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_region);
patch([x fliplr(x)], [y1 fliplr(y2)],color_region,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_region,'FaceAlpha',0.2,'EdgeColor','none')
hold on
y1=surr_s(1,:);
y2=surr_s(2,:);
y3=surr_s(3,:);
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_surr);
patch([x fliplr(x)], [y1 fliplr(y2)],color_surr,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_surr,'FaceAlpha',0.2,'EdgeColor','none')
ax2=subplot(2,2,3)
y1=med_thal_s{2}(1,:);
y2=med_thal_s{2}(2,:);
y3=med_thal_s{2}(3,:);
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_region);
patch([x fliplr(x)], [y1 fliplr(y2)],color_region,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_region,'FaceAlpha',0.2,'EdgeColor','none')
hold on
y1=surr_s(1,:);
y2=surr_s(2,:);
y3=surr_s(3,:);
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_surr);
patch([x fliplr(x)], [y1 fliplr(y2)],color_surr,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_surr,'FaceAlpha',0.2,'EdgeColor','none')


%-------SNR

clearvars -except ax1 ax2
cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A1_Thal/mat')
load ('SNR_fig','med_contacts','surr','surr_s')
out3_real= med_contacts;
clear med_contacts
out3_sur= surr;
clear surr
color_BZ=[0.5 0 0.5];
color_SNR=[0 0 0.5];
color_surr=[0.3 0.3 0.3];
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
    ax3=subplot(2,2,2)
    patch([beg(1) beg(end) beg(end) beg(1)]-200,[0.2 0.2 0.225 0.225],color_region,'EdgeColor','none')
end
hold on
clear beg
beg=find(st(2,:)<0.05 & st2(2,:)~=0);
if ~isempty(beg)
    beg(1)
    beg(find(diff(beg)>1))
    beg(find(diff(beg)>1)+1)
    beg(end)
    ax4=subplot(2,2,4)
    patch([beg(1) beg(end) beg(end) beg(1)]-200,[0.2 0.2 0.225 0.225],color_region,'EdgeColor','none')
end
hold on
load ('SNR_fig','med_thal_s')
ax3=subplot(2,2,2)
x=[-500:500];
y1=med_thal_s{1}(1,:);
y2=med_thal_s{1}(2,:);
y3=med_thal_s{1}(3,:);
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_region);
patch([x fliplr(x)], [y1 fliplr(y2)],color_region,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_region,'FaceAlpha',0.2,'EdgeColor','none')
hold on
y1=surr_s(1,:);
y2=surr_s(2,:);
y3=surr_s(3,:);
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_surr);
patch([x fliplr(x)], [y1 fliplr(y2)],color_surr,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_surr,'FaceAlpha',0.2,'EdgeColor','none')
ax4=subplot(2,2,4)
y1=med_thal_s{2}(1,:);
y2=med_thal_s{2}(2,:);
y3=med_thal_s{2}(3,:);
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_region);
patch([x fliplr(x)], [y1 fliplr(y2)],color_region,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_region,'FaceAlpha',0.2,'EdgeColor','none')
hold on
y1=surr_s(1,:);
y2=surr_s(2,:);
y3=surr_s(3,:);
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_surr);
patch([x fliplr(x)], [y1 fliplr(y2)],color_surr,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_surr,'FaceAlpha',0.2,'EdgeColor','none')
ylim([ax1 ax2 ax3 ax4],[-0.25 0.25])
xticks([ax1 ax2 ax3 ax4],[-500:250:500])
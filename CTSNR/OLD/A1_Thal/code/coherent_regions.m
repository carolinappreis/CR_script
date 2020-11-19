%Ploting median change in beta activity at the subcortical level aligned with short bursts (only beta coherent)
%significant different between the change in beta evoked by bursts vs.
%change in beta evoked by surrogates (function hayriye_c needed).

clear all
cd ('\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A1_Thal\mat')
load ('BZ_fig_new','all_contacts','surr2','surr_s')
out3_real= all_contacts;
clear all_contacts
out3_sur= surr2;
clear surr2
color_BZ=[0.5 0 0.5];
color_SNR=[0 0 0.5];
color_surr=[0.3 0.3 0.3];
color_region=color_BZ;


st=NaN(2,551);
for ind=1:2; clear A; A(1:size(out3_real,2),1:551)=out3_real(ind,:,300:850);
    clear B; B(1:size(out3_sur,1),1:551)=out3_sur(:,300:850);
    hayriye_c; st(ind,:)=stats.prob; st2(ind,:)=stats.posclusterslabelmat; % needs to seth path to C:\Users\creis\Documents\MATLAB\spm12_clean_noconf\external\fieldtrip
%     figure(ind)
%     plot(median(B))
%     hold on
%     plot(median(A))
%     plot(st(ind,:),'r.')
%     plot(st2(ind,:),'b.')
end
ax1=subplot(2,1,1)
title ('BZ','FontSize',14)
beg=find(st(1,:)<0.05 & st2(1,:)~=0);
if ~isempty(beg)
    beg(1)
    beg(find(diff(beg)>1))
    beg(find(diff(beg)>1)+1)
    beg(end)
    figure(1)
    patch([beg(1) beg(end) beg(end) beg(1)]-200,[0.15 0.15 0.175 0.175],color_region+[0 0 0.3],'EdgeColor','none')
end
hold on
clear beg
beg=find(st(2,:)<0.05 & st2(2,:)~=0);
if ~isempty(beg)
    beg(1)
    beg(find(diff(beg)>1))
    beg(find(diff(beg)>1)+1)
    beg(end)
    patch([beg(1) beg(end) beg(end) beg(1)]-200,[0.180 0.180 0.205 0.205],color_region - [0.3 0 0.3],'EdgeColor','none')
end
hold on
load ('BZ_fig_new','med_thal_s')
x=[-500:500];
y1=med_thal_s{1}(1,:);
y2=med_thal_s{1}(2,:);
y3=med_thal_s{1}(3,:);
plot(x,y2,'DisplayName','Beta change with short bursts')
set(plot(x,y2),'LineWidth',1.5,'Color',color_region+[0 0 0.3]);
patch([x fliplr(x)], [y1 fliplr(y2)],color_region+[0 0 0.3],'FaceAlpha',0.1,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_region+[0 0 0.3],'FaceAlpha',0.1,'EdgeColor','none')
hold on
y1=surr_s(1,:);
y2=surr_s(2,:);
y3=surr_s(3,:);
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_surr);
patch([x fliplr(x)], [y1 fliplr(y2)],color_surr,'FaceAlpha',0.1,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_surr,'FaceAlpha',0.1,'EdgeColor','none')
hold on
y1=med_thal_s{2}(1,:);
y2=med_thal_s{2}(2,:);
y3=med_thal_s{2}(3,:);
plot(x,y2,'DisplayName','Beta change with long bursts')
set(plot(x,y2),'LineWidth',1.5,'Color',color_region - [0.3 0 0.3]);
patch([x fliplr(x)], [y1 fliplr(y2)],color_region- [0.3 0 0.3],'FaceAlpha',0.1,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_region- [0.3 0 0.3],'FaceAlpha',0.1,'EdgeColor','none')



%%%-------SNR


load ('SNR_fig_new','all_contacts','surr2','surr_s')
out3_real= all_contacts;
clear all_contacts
out3_sur= surr2;
clear surr2
color_BZ=[0.5 0 0.5];
color_SNR=[0 0 0.5];
color_surr=[0.3 0.3 0.3];
color_region=color_SNR;

st=NaN(2,551);
for ind=1:2; clear A; A(1:size(out3_real,2),1:551)=out3_real(ind,:,300:850);
    clear B; B(1:size(out3_sur,1),1:551)=out3_sur(:,300:850);
  hayriye_c; st(ind,:)=stats.prob; st2(ind,:)=stats.posclusterslabelmat;
%     figure(ind)
%     plot(median(B))
%     hold on
%     plot(median(A))
%     plot(st(ind,:),'r.')
%     plot(st2(ind,:),'b.')
end

ax2=subplot(2,1,2)
title ('SNr','FontSize',14)
beg=find(st(1,:)<0.05 & st2(1,:)~=0);
if ~isempty(beg)
    beg(1)
    beg(find(diff(beg)>1))
    beg(find(diff(beg)>1)+1)
    beg(end)
    patch([beg(1) beg(end) beg(end) beg(1)]-200,[0.15 0.15 0.175 0.175],color_region+[0 0 0.3],'EdgeColor','none')
end
hold on
clear beg
beg=find(st(2,:)<0.05 & st2(2,:)~=0);
if ~isempty(beg)
    beg(1)
    beg(find(diff(beg)>1))
    beg(find(diff(beg)>1)+1)
    beg(end)
    patch([beg(1) beg(end) beg(end) beg(1)]-200,[0.180 0.180 0.205 0.205],color_region- [0 0 0.3],'EdgeColor','none')
end
hold on
load ('SNR_fig_new','med_thal_s')
x=[-500:500];
y1=med_thal_s{1}(1,:);
y2=med_thal_s{1}(2,:);
y3=med_thal_s{1}(3,:);
plot(x,y2,'DisplayName','Beta change with short bursts')
set(plot(x,y2),'LineWidth',1.5,'Color',color_region+[0 0 0.3]);
patch([x fliplr(x)], [y1 fliplr(y2)],color_region+[0 0 0.3],'FaceAlpha',0.1,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_region+[0 0 0.3],'FaceAlpha',0.1,'EdgeColor','none')
hold on
y1=surr_s(1,:);
y2=surr_s(2,:);
y3=surr_s(3,:);
plot(x,y2,'DisplayName','Beta change outside bursts')
set(plot(x,y2),'LineWidth',1.5,'Color',color_surr);
patch([x fliplr(x)], [y1 fliplr(y2)],color_surr,'FaceAlpha',0.1,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_surr,'FaceAlpha',0.1,'EdgeColor','none')
y1=med_thal_s{2}(1,:);
y2=med_thal_s{2}(2,:);
y3=med_thal_s{2}(3,:);
plot(x,y2,'DisplayName','Beta change with long bursts')
set(plot(x,y2),'LineWidth',1.5,'Color',color_region- [0 0 0.2]);
patch([x fliplr(x)], [y1 fliplr(y2)],color_region- [0 0 0.3],'FaceAlpha',0.1,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_region- [0 0 0.3],'FaceAlpha',0.1,'EdgeColor','none')
xlabel( 'msec','FontSize',14)
ylabel( 'Beta amplitude (%Baseline)','FontSize',14)

ylim( [ax1 ax2] , [-0.22 0.22])
yticks ( [ax1 ax2] , [-0.2 -0.1 0 0.1 0.2 ])
yticklabels ( [ax1 ax2] , {'-20' ,'-10' ,'0' ,'10' ,'20'})




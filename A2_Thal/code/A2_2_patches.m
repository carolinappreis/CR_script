
clear all
% cd ('\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A2_Thal\mat')
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A2_Thal/mat')
load ('A2_2_BZ.mat')
figure(1)
subplot(2,1,1)
x=[-200:200];
y1=zscore_best_stats(1,select);
y2=zscore_best_stats(2,select);
y3=zscore_best_stats(3,select);
load ('A2_2_SNR.mat')
y4=zscore_best_stats(1,select);
y5=zscore_best_stats(2,select);
y6=zscore_best_stats(3,select);
color_bz=[0.5 0 0.5];
color_snr=[0 0 0.5];
plot(x, y2, 'DisplayName','BZ aligned to short bursts')
set(plot(x,y2),'LineWidth',1.5,'Color',color_bz);
patch([x fliplr(x)], [y1 fliplr(y2)],[color_bz],'FaceAlpha',[0.2],'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],[color_bz],'FaceAlpha',[0.2],'EdgeColor','none')
hold on
plot(x, y5)
set(plot(x,y5),'LineWidth',1.5,'Color',color_snr);
patch([x fliplr(x)], [y4 fliplr(y5)],[color_snr],'FaceAlpha',[0.2],'EdgeColor','none')
patch([x fliplr(x)], [y5 fliplr(y6)],[color_snr],'FaceAlpha',[0.2],'EdgeColor','none')
title ('Z-score of oscillatory response phase-locked to cortical beta bursts','FontSize',12)
xlabel ('msec','FontSize',12)
ylabel ('zscore','FontSize',12)
ylim ([-5 5])

subplot(2,1,2)
clear all
cd ('\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A2_Thal\mat')
load ('A2_2_BZ_surr.mat')
x=[-200:200];
y1=zscore_best_stats(1,select);
y2=zscore_best_stats(2,select);
y3=zscore_best_stats(3,select);
load ('A2_2_SNR_surr.mat')
y4=zscore_best_stats(1,select);
y5=zscore_best_stats(2,select);
y6=zscore_best_stats(3,select);
color_bz=[0.5 0 0.5];
color_snr=[0 0 0.5];
plot(x, y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_bz);
patch([x fliplr(x)], [y1 fliplr(y2)],[color_bz],'FaceAlpha',[0.2],'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],[color_bz],'FaceAlpha',[0.2],'EdgeColor','none')
hold on
plot(x, y5)
set(plot(x,y5),'LineWidth',1.5,'Color',color_snr);
patch([x fliplr(x)], [y4 fliplr(y5)],[color_snr],'FaceAlpha',[0.2],'EdgeColor','none')
patch([x fliplr(x)], [y5 fliplr(y6)],[color_snr],'FaceAlpha',[0.2],'EdgeColor','none')
title ('Z-score of oscillatory response non phase-locked to cortical beta bursts','FontSize',12)
xlabel ('msec','FontSize',12)
ylabel ('zscore','FontSize',12)
ylim ([-5 5])

clear all
cd ('\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A2_Thal\mat')
figure(2)
subplot(2,1,1)
load ('A2_2_BZ.mat')
x=[-200:200];
y1=zbest_stats(1,:);
y2=zbest_stats(2,:);
y3=zbest_stats(3,:);
color_bz=[0.5 0 0.5];
plot(x, y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_bz);
patch([x fliplr(x)], [y1 fliplr(y2)],[color_bz],'FaceAlpha',[0.2],'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],[color_bz],'FaceAlpha',[0.2],'EdgeColor','none')
hold on
clear all
load ('A2_2_SNR.mat')
x=[-200:200];
y1=zbest_stats(1,:);
y2=zbest_stats(2,:);
y3=zbest_stats(3,:);
color_snr=[0 0 0.5];
plot(x, y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_snr);
patch([x fliplr(x)], [y1 fliplr(y2)],[color_snr],'FaceAlpha',[0.2],'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],[color_snr],'FaceAlpha',[0.2],'EdgeColor','none')
title ('Envelope of z-scored responses phase-locked to ctx beta bursts','FontSize',12)
xlabel ('msec','FontSize',12)
ylabel ('envelope zscore beta evoked','FontSize',12)
ylim ([0 6])
clear all

subplot(2,1,2)
load ('A2_2_BZ_surr.mat')
x=[-200:200];
y1=zbest_stats(1,:);
y2=zbest_stats(2,:);
y3=zbest_stats(3,:);
color_bz=[0.5 0 0.5];
plot(x, y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_bz);
patch([x fliplr(x)], [y1 fliplr(y2)],[color_bz],'FaceAlpha',[0.2],'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],[color_bz],'FaceAlpha',[0.2],'EdgeColor','none')
hold on
clear all
load ('A2_2_SNR_surr.mat')
x=[-200:200];
y1=zbest_stats(1,:);
y2=zbest_stats(2,:);
y3=zbest_stats(3,:);
color_snr=[0 0 0.5];
plot(x, y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_snr);
patch([x fliplr(x)], [y1 fliplr(y2)],[color_snr],'FaceAlpha',[0.2],'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],[color_snr],'FaceAlpha',[0.2],'EdgeColor','none')
title ('Envelope of z-scored responses non phase-locked to ctx beta bursts','FontSize',12)
xlabel ('msec','FontSize',12)
ylabel ('envelope zscore beta evoked','FontSize',12)
ylim ([0 6])



% clear all
% cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A2_Thal/mat')
% load ('A2_2_peaks_all.mat')
% figure(3)
% x=[peak_bz_snr(1,1):1:peak_bz_snr(1,3)];
% y1=3.5*ones(length(peak_bz_snr(1,1):peak_bz_snr(1,3)),1)';
% y2=3.3*ones(length(peak_bz_snr(1,1):peak_bz_snr(1,3)),1)';
% patch( [x fliplr(x)],[y1 fliplr(y2)],[0.5 0 0.5],'FaceAlpha',[0.2],'EdgeColor','none')
% hold on
% plot ([peak_bz_snr(1,2) peak_bz_snr(1,2)],[3.4 3.4],'o',...
%     'MarkerFaceColor',[0.5 0 0.5],'MarkerEdgeColor','none','MarkerSize',5)
% hold on
% a=[peak_bz_snr(2,1):1:peak_bz_snr(2,3)];
% b=2.5*ones(length(peak_bz_snr(2,1):peak_bz_snr(2,3)),1)';
% c=2.3*ones(length(peak_bz_snr(2,1):peak_bz_snr(2,3)),1)';
% patch( [a fliplr(a)],[b fliplr(c)],[0 0 0.5],'FaceAlpha',[0.2],'EdgeColor','none')
% hold on
% plot ([peak_bz_snr(2,2) peak_bz_snr(2,2)],[2.4 2.4],'o',...
%     'MarkerFaceColor',[0 0 0.5],'MarkerEdgeColor','none','MarkerSize',5)
% ylim ([1.5 4])
% xlim ([0 300])
% yticks ([2.4 3.4])
% yticklabels ({'SNr','BZ'})
% xlabel ('msec','FontSize',12)

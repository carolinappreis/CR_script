
clear all
cd C:\Users\wolf5173\Documents\Carolina_code\codes_thal\A2_Thal\mat\A2
load ('A2_evoked_fig.mat')
figure(1)
subplot(2,1,1)
x=[-200:200];
y1=filt_best_stats_bz(1,select);
y2=filt_best_stats_bz(2,select);
y3=filt_best_stats_bz(3,select);
y4=filt_best_stats_snr(1,select);
y5=filt_best_stats_snr(2,select);
y6=filt_best_stats_snr(3,select);
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
title ('Phase locked thalamic response to cortical beta burst ','FontSize',12)
xlabel ('msec','FontSize',10)
ylabel ('Amplitude(mV)','FontSize',10)

subplot(2,1,2)
clear all
cd C:\Users\wolf5173\Documents\Carolina_code\codes_thal\A2_Thal\mat\A2
load ('A2_evoked_surr_fig.mat')
x=[-200:200];
y1=filt_best_stats_bz(1,select);
y2=filt_best_stats_bz(2,select);
y3=filt_best_stats_bz(3,select);
y4=filt_best_stats_snr(1,select);
y5=filt_best_stats_snr(2,select);
y6=filt_best_stats_snr(3,select);
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
title ('Phase locked thalamic response to surrogates ','FontSize',12)
xlabel ('msec','FontSize',10)
ylabel ('Amplitude(mV)','FontSize',10)

clear all
cd C:\Users\wolf5173\Documents\Carolina_code\codes_thal\A2_Thal\mat\A2
figure(2)
subplot(2,1,1)
load ('A2_BZ.mat')
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
load ('A2_SNR.mat')
x=[-200:200];
y1=zbest_stats(1,:);
y2=zbest_stats(2,:);
y3=zbest_stats(3,:);
color_snr=[0 0 0.5];
plot(x, y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_snr);
patch([x fliplr(x)], [y1 fliplr(y2)],[color_snr],'FaceAlpha',[0.2],'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],[color_snr],'FaceAlpha',[0.2],'EdgeColor','none')
title ('Z-score of beta thalamic evelope locked to ctx beta bursts','FontSize',12)
xlabel ('msec','FontSize',10)
ylabel ('Amplitude(mV)','FontSize',10)
clear all

subplot(2,1,2)
load ('A2_BZ_surr.mat')
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
load ('A2_SNR_surr.mat')
x=[-200:200];
y1=zbest_stats(1,:);
y2=zbest_stats(2,:);
y3=zbest_stats(3,:);
color_snr=[0 0 0.5];
plot(x, y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_snr);
patch([x fliplr(x)], [y1 fliplr(y2)],[color_snr],'FaceAlpha',[0.2],'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],[color_snr],'FaceAlpha',[0.2],'EdgeColor','none')
title ('Z-score envelope of phase locked thalamic response to surrogates','FontSize',12)
xlabel ('msec','FontSize',10)
ylabel ('Z-score beta envelope','FontSize',10)


clear all
cd C:\Users\wolf5173\Documents\Carolina_code\codes_thal\A2_Thal\mat
load ('peaks_all.mat')
figure(3)
x=[peak_bz_snr(1,1):1:peak_bz_snr(1,3)];
y1=3.5*ones(length(peak_bz_snr(1,1):peak_bz_snr(1,3)),1)';
y2=3.3*ones(length(peak_bz_snr(1,1):peak_bz_snr(1,3)),1)';
patch( [x fliplr(x)],[y1 fliplr(y2)],[0.5 0 0.5],'FaceAlpha',[0.2],'EdgeColor','none')
hold on
plot ([peak_bz_snr(1,2) peak_bz_snr(1,2)],[3.4 3.4],'o',...
    'MarkerFaceColor',[0.5 0 0.5],'MarkerEdgeColor','none','MarkerSize',5)
hold on
a=[peak_bz_snr(2,1):1:peak_bz_snr(2,3)];
b=2.5*ones(length(peak_bz_snr(2,1):peak_bz_snr(2,3)),1)';
c=2.3*ones(length(peak_bz_snr(2,1):peak_bz_snr(2,3)),1)';
patch( [a fliplr(a)],[b fliplr(c)],[0 0 0.5],'FaceAlpha',[0.2],'EdgeColor','none')
hold on
plot ([peak_bz_snr(2,2) peak_bz_snr(2,2)],[2.4 2.4],'o',...
    'MarkerFaceColor',[0 0 0.5],'MarkerEdgeColor','none','MarkerSize',5)
ylim ([1.5 4])
xlim ([0 300])
yticks ([2.4 3.4])
yticklabels ({'SNr','BZ'})
xlabel ('msec','FontSize',10)

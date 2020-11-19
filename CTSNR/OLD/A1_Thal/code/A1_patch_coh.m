clear all
cd C:\Users\wolf5173\Documents\Carolina_code\codes_thal\A1_Thal\mat

%50-median ctx bursts
load 'BZ.mat'
figure(1)
subplot(3,2,1)
x=[-500:500];
y1=ctx_sl{1}(1,:);
y2=ctx_sl{1}(2,:);
y3=ctx_sl{1}(3,:);
color_bz=[0.5 0 0.5];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_bz);
patch([x fliplr(x)], [y1 fliplr(y2)],color_bz,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_bz,'FaceAlpha',0.2,'EdgeColor','none')
hold on
clear all
load 'SNR.mat'
x=[-500:500];
y1=ctx_sl{1}(1,:);
y2=ctx_sl{1}(2,:);
y3=ctx_sl{1}(3,:);
color_snr=[0 0 0.5];
plot(x, y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_snr);
patch([x fliplr(x)], [y1 fliplr(y2)],[color_snr],'FaceAlpha',[0.2],'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],[color_snr],'FaceAlpha',[0.2],'EdgeColor','none')
title ('short cortical burst','FontSize',12)
xlabel ('msec','FontSize',10)
ylabel ('Power(mV^2)','FontSize',10)
xticks([-500:250:500])



%200p ctx bursts
load 'BZ.mat'
subplot(3,2,2)
x=[-500:500];
y1=ctx_sl{2}(1,:);
y2=ctx_sl{2}(2,:);
y3=ctx_sl{2}(3,:);
color_bz=[0.5 0 0.5];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_bz);
patch([x fliplr(x)], [y1 fliplr(y2)],color_bz,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_bz,'FaceAlpha',0.2,'EdgeColor','none')
hold on
clear all
load 'SNR.mat'
x=[-500:500];
y1=ctx_sl{2}(1,:);
y2=ctx_sl{2}(2,:);
y3=ctx_sl{2}(3,:);
color_snr=[0 0 0.5];
plot(x, y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_snr);
patch([x fliplr(x)], [y1 fliplr(y2)],[color_snr],'FaceAlpha',[0.2],'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],[color_snr],'FaceAlpha',[0.2],'EdgeColor','none')
title ('long cortical burst','FontSize',12)
xlabel ('msec','FontSize',10)
ylabel ('Power(mV^2)','FontSize',10)
xticks([-500:250:500])

%more than median duration bursts
load 'BZ.mat'
subplot(3,2,3)
x=[-500:500];
y1=thal_sls{1}(1,:);
y2=thal_sls{1}(2,:);
y3=thal_sls{1}(3,:);
color_bz=[0.5 0 0.5];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_bz);
patch([x fliplr(x)], [y1 fliplr(y2)],color_bz,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_bz,'FaceAlpha',0.2,'EdgeColor','none')
hold on
clear all
load 'SNR.mat'
x=[-500:500];
y1=thal_sls{1}(1,:);
y2=thal_sls{1}(2,:);
y3=thal_sls{1}(3,:);
color_snr=[0 0 0.5];
plot(x, y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_snr);
patch([x fliplr(x)], [y1 fliplr(y2)],[color_snr],'FaceAlpha',[0.2],'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],[color_snr],'FaceAlpha',[0.2],'EdgeColor','none')
title ('thalamic beta spectral change','FontSize',12)
xlabel ('msec','FontSize',10)
ylabel ('Power(mV^2)','FontSize',10)
xticks([-500:250:500])


%200plus bursts
subplot(3,2,4)
load 'BZ.mat'
x=[-500:500];
y1=thal_sls{2}(1,:);
y2=thal_sls{2}(2,:);
y3=thal_sls{2}(3,:);
color_bz=[0.5 0 0.5];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_bz);
patch([x fliplr(x)], [y1 fliplr(y2)],color_bz,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_bz,'FaceAlpha',0.2,'EdgeColor','none')
hold on
clear all
load 'SNR.mat'
x=[-500:500];
y1=thal_sls{2}(1,:);
y2=thal_sls{2}(2,:);
y3=thal_sls{2}(3,:);
color_snr=[0 0 0.5];
plot(x, y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_snr);
patch([x fliplr(x)], [y1 fliplr(y2)],[color_snr],'FaceAlpha',[0.2],'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],[color_snr],'FaceAlpha',[0.2],'EdgeColor','none')
title ('thalamic beta spectral change','FontSize',12)
xlabel ('msec','FontSize',10)
ylabel ('Power(mV^2)','FontSize',10)
xticks([-500:250:500])


%surrogates
load 'BZ.mat'
subplot(3,2,[5 6])
x=[-500:500];
y1=thal_sls{3}(1,:);
y2=thal_sls{3}(2,:);
y3=thal_sls{3}(3,:);
color_bz=[0.5 0 0.5];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_bz);
patch([x fliplr(x)], [y1 fliplr(y2)],color_bz,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_bz,'FaceAlpha',0.2,'EdgeColor','none')
hold on
clear all
load 'SNR.mat'
x=[-500:500];
y1=thal_sls{3}(1,:);
y2=thal_sls{3}(2,:);
y3=thal_sls{3}(3,:);
color_snr=[0 0 0.5];
plot(x, y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_snr);
patch([x fliplr(x)], [y1 fliplr(y2)],[color_snr],'FaceAlpha',[0.2],'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],[color_snr],'FaceAlpha',[0.2],'EdgeColor','none')
title ('Surrogates','FontSize',12)
xlabel ('msec','FontSize',10)
ylabel ('Power(mV^2)','FontSize',10)
xticks([-500:100:500])
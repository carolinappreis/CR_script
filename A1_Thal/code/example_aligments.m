% change BZ.mat vs. SNR.mat and color_BZ vs. color_SNR
clear all
cd C:\Users\wolf5173\Documents\Carolina_code\codes_thal\A1_Thal\mat
figure(1)
subplot(2,2,1)
load 'SNR.mat'
x=[-500:500];
x=[-500:500];
y1=ctx_sl{1}(1,:);
y2=ctx_sl{1}(2,:);
y3=ctx_sl{1}(3,:);
color_SNR=[0 0 0.5];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_SNR);
patch([x fliplr(x)], [y1 fliplr(y2)],color_SNR,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_SNR,'FaceAlpha',0.2,'EdgeColor','none')
xticks([-500:250:500])
title ('short cortical burst','FontSize',12)
xlabel ('msec','FontSize',10)
ylabel ('Power(mV^2)','FontSize',10)

subplot(2,2,3)
load 'SNR.mat'
x=[-500:500];
y1=thal_sls{1}(1,:);
y2=thal_sls{1}(2,:);
y3=thal_sls{1}(3,:);
color_SNR=[0 0 0.5];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_SNR);
patch([x fliplr(x)], [y1 fliplr(y2)],color_SNR,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_SNR,'FaceAlpha',0.2,'EdgeColor','none')
xticks([-500:250:500])
title ('SNR power changes','FontSize',12)
xlabel ('msec','FontSize',10)
ylabel ('Power(mV^2)','FontSize',10)

subplot(2,2,2)
load 'SNR.mat'
x=[-500:500];
y1=ctx_sl{2}(1,:);
y2=ctx_sl{2}(2,:);
y3=ctx_sl{2}(3,:);
color_SNR=[0 0 0.5];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_SNR);
patch([x fliplr(x)], [y1 fliplr(y2)],color_SNR,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_SNR,'FaceAlpha',0.2,'EdgeColor','none')
xticks([-500:250:500])
title ('long cortical burst','FontSize',12)
xlabel ('msec','FontSize',10)
ylabel ('Power(mV^2)','FontSize',10)

subplot(2,2,4)
load 'SNR.mat'
x=[-500:500];
y1=thal_sls{2}(1,:);
y2=thal_sls{2}(2,:);
y3=thal_sls{2}(3,:);
color_SNR=[0 0 0.5];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_SNR);
patch([x fliplr(x)], [y1 fliplr(y2)],color_SNR,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_SNR,'FaceAlpha',0.2,'EdgeColor','none')
xticks([-500:250:500])
title ('SNR power changes','FontSize',12)
xlabel ('msec','FontSize',10)
ylabel ('Power(mV^2)','FontSize',10)
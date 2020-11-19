clear all
cd C:\Users\wolf5173\Documents\Carolina_code\codes_thal\A1_Thal\mat

figure('InvertHardcopy','off','Color',[1 1 1],...
    'OuterPosition',[433 519 581.333333333333 352]);

load 'nc_BZ_50_200.mat'
x=[-500:500];
y1=nc_BZ_50_200_thal(1,:);
y2=nc_BZ_50_200_thal(2,:);
y3=nc_BZ_50_200_thal(3,:);
color_BZ=[0.5 0 0.5];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_BZ);
patch([x fliplr(x)], [y1 fliplr(y2)],color_BZ,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_BZ,'FaceAlpha',0.2,'EdgeColor','none')
hold on

clear all
load 'nc_SNR_50_200.mat'
x=[-500:500];
y1=nc_SNR_50_200_thal(1,:);
y2=nc_SNR_50_200_thal(2,:);
y3=nc_SNR_50_200_thal(3,:);
color_SNR=[0 0 0.5];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_SNR);
patch([x fliplr(x)], [y1 fliplr(y2)],color_SNR,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_SNR,'FaceAlpha',0.2,'EdgeColor','none')

clear all
load 'nc_ZI_50_200.mat'
x=[-500:500];
y1=nc_ZI_50_200_thal(1,:);
y2=nc_ZI_50_200_thal(2,:);
y3=nc_ZI_50_200_thal(3,:);
color_ZI=[0.5 0 0];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_ZI);
patch([x fliplr(x)], [y1 fliplr(y2)],color_ZI,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_ZI,'FaceAlpha',0.2,'EdgeColor','none')

clear all
load 'nc_AV_50_200.mat'
x=[-500:500];
y1=nc_AV_50_200_thal(1,:);
y2=nc_AV_50_200_thal(2,:);
y3=nc_AV_50_200_thal(3,:);
color_AV=[0 0.5 0];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_AV);
patch([x fliplr(x)], [y1 fliplr(y2)],color_AV,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_AV,'FaceAlpha',0.2,'EdgeColor','none')


% clear all
% load 'nc_AM_50_200.mat'
% x=[-500:500];
% y1=nc_AM_50_200_thal(1,:);
% y2=nc_AM_50_200_thal(2,:);
% y3=nc_AM_50_200_thal(3,:);
% color_AM=[0.5 0 0.5];
% plot(x,y2)
% set(plot(x,y2),'LineWidth',1.5,'Color',color_AM);
% patch([x fliplr(x)], [y1 fliplr(y2)],color_AM,'FaceAlpha',0.2,'EdgeColor','none')
% patch([x fliplr(x)], [y2 fliplr(y3)],color_AM,'FaceAlpha',0.2,'EdgeColor','none')
% 

clear all
load 'nc_RT_50_200.mat'
x=[-500:500];
y1=nc_RT_50_200_thal(1,:);
y2=nc_RT_50_200_thal(2,:);
y3=nc_RT_50_200_thal(3,:);
color_RT=[0.5 0.5 0.5];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_RT);
patch([x fliplr(x)], [y1 fliplr(y2)],color_RT,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_RT,'FaceAlpha',0.2,'EdgeColor','none')
hold on

clear all
load 'nc_NW_50_200.mat'
x=[-500:500];
y1=nc_NW_50_200_thal(1,:);
y2=nc_NW_50_200_thal(2,:);
y3=nc_NW_50_200_thal(3,:);
color_NW=[0 0.1 0.2];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_NW);
patch([x fliplr(x)], [y1 fliplr(y2)],color_NW,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_NW,'FaceAlpha',0.2,'EdgeColor','none')

clear all
load 'nc_PC_50_200.mat'
x=[-500:500];
y1=nc_PC_50_200_thal(1,:);
y2=nc_PC_50_200_thal(2,:);
y3=nc_PC_50_200_thal(3,:);
color_PC=[0.4 0 0];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_PC);
patch([x fliplr(x)], [y1 fliplr(y2)],color_PC,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_PC,'FaceAlpha',0.2,'EdgeColor','none')

hold on
clear all
load 'nc_CZ_50_200.mat'
x=[-500:500];
y1=nc_CZ_50_200_thal(1,:);
y2=nc_CZ_50_200_thal(2,:);
y3=nc_CZ_50_200_thal(3,:);
color_CZ=[0 0.2 0.5];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_CZ);
patch([x fliplr(x)], [y1 fliplr(y2)],color_CZ,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_CZ,'FaceAlpha',0.2,'EdgeColor','none')

title ('Thalamic beta spectral changes aligned to ctx short beta bursts','FontSize',12)
xlabel('msec','FontSize',12)
ylabel('Power(mV^2)','FontSize',12)
clear all
cd C:\Users\wolf5173\Documents\Carolina_code\codes_thal\A1_Thal\mat

load 'nc_BZ_200p.mat'
figure1 = figure('InvertHardcopy','off','Color',[1 1 1],...
    'OuterPosition',[1009 27.6666666666667 373.333333333333 843.333333333333]);
subplot(8,1,1)
x=[-500:500];
y1=nc_BZ_200p_thal(1,:);
y2=nc_BZ_200p_thal(2,:);
y3=nc_BZ_200p_thal(3,:);
color_BZ=[0.5 0 0.5];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_BZ);
patch([x fliplr(x)], [y1 fliplr(y2)],color_BZ,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_BZ,'FaceAlpha',0.2,'EdgeColor','none')
title ('BZ')
hold on

clear all
load 'nc_SNR_200p.mat'
subplot(8,1,2)
x=[-500:500];
y1=nc_SNR_200p_thal(1,:);
y2=nc_SNR_200p_thal(2,:);
y3=nc_SNR_200p_thal(3,:);
color_SNR=[0 0 0.5];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_SNR);
patch([x fliplr(x)], [y1 fliplr(y2)],color_SNR,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_SNR,'FaceAlpha',0.2,'EdgeColor','none')
title ('SNr')

clear all
load 'nc_CZ_200p.mat'
subplot(8,1,3)
x=[-500:500];
y1=nc_CZ_200p_thal(1,:);
y2=nc_CZ_200p_thal(2,:);
y3=nc_CZ_200p_thal(3,:);
color_CZ=[0 0.2 0.5];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_CZ);
patch([x fliplr(x)], [y1 fliplr(y2)],color_CZ,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_CZ,'FaceAlpha',0.2,'EdgeColor','none')
title ('CZ')


clear all
load 'nc_AV_200p.mat'
subplot(8,1,4)
x=[-500:500];
y1=nc_AV_200p_thal(1,:);
y2=nc_AV_200p_thal(2,:);
y3=nc_AV_200p_thal(3,:);
color_AV=[0 0.5 0];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_AV);
patch([x fliplr(x)], [y1 fliplr(y2)],color_AV,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_AV,'FaceAlpha',0.2,'EdgeColor','none')
title ('AV')

% clear all
% load 'nc_AM_200p.mat'
% x=[-500:500];
% y1=nc_AM_200p_thal(1,:);
% y2=nc_AM_200p_thal(2,:);
% y3=nc_AM_200p_thal(3,:);
% color_AM=[0.5 0 0.5];
% plot(x,y2)
% set(plot(x,y2),'LineWidth',1.5,'Color',color_AM);
% patch([x fliplr(x)], [y1 fliplr(y2)],color_AM,'FaceAlpha',0.2,'EdgeColor','none')
% patch([x fliplr(x)], [y2 fliplr(y3)],color_AM,'FaceAlpha',0.2,'EdgeColor','none')
% 

clear all
load 'nc_RT_200p.mat'
subplot(8,1,5)
x=[-500:500];
y1=nc_RT_200p_thal(1,:);
y2=nc_RT_200p_thal(2,:);
y3=nc_RT_200p_thal(3,:);
color_RT=[0.5 0.5 0.5];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_RT);
patch([x fliplr(x)], [y1 fliplr(y2)],color_RT,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_RT,'FaceAlpha',0.2,'EdgeColor','none')
title ('Rt')


clear all
load 'nc_NW_200p.mat'
subplot(8,1,6)
x=[-500:500];
y1=nc_NW_200p_thal(1,:);
y2=nc_NW_200p_thal(2,:);
y3=nc_NW_200p_thal(3,:);
color_NW=[0 0.1 0.2];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_NW);
patch([x fliplr(x)], [y1 fliplr(y2)],color_NW,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_NW,'FaceAlpha',0.2,'EdgeColor','none')
title ('NW')


clear all
load 'nc_PC_200p.mat'
subplot(8,1,7)
x=[-500:500];
y1=nc_PC_200p_thal(1,:);
y2=nc_PC_200p_thal(2,:);
y3=nc_PC_200p_thal(3,:);
color_PC=[0.4 0 0];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_PC);
patch([x fliplr(x)], [y1 fliplr(y2)],color_PC,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_PC,'FaceAlpha',0.2,'EdgeColor','none')
title ('PC')

clear all
load 'nc_ZI_200p.mat'
subplot(8,1,8)
x=[-500:500];
y1=nc_ZI_200p_thal(1,:);
y2=nc_ZI_200p_thal(2,:);
y3=nc_ZI_200p_thal(3,:);
color_ZI=[0.5 0 0];
plot(x,y2)
set(plot(x,y2),'LineWidth',1.5,'Color',color_ZI);
patch([x fliplr(x)], [y1 fliplr(y2)],color_ZI,'FaceAlpha',0.2,'EdgeColor','none')
patch([x fliplr(x)], [y2 fliplr(y3)],color_ZI,'FaceAlpha',0.2,'EdgeColor','none')
title ('ZI')
xlabel('msec','FontSize',12)


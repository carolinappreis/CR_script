clear all
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal')
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal')

load('BZ_opt.mat');  
load('SNR_opt.mat'); 
color_s=[0 0 0];
color_ctxb=[0.5 0.5 0.5];
time=1:501;


for m=1:size(SNR.idrat,1)
    r=SNR.idrat(m);
new{m,1}=(SNR.power_thal{r,:});
end

fig=figure();
i=vertcat(new{:}); color_b= [0 0 0.5];
% i=SNR.power_ctx; color_b=[0.5 0.5 0.5];
 
y2=log(mean(i,1));
y1=log(mean(i,1)+std(i)./sqrt(size(i,1)));
y3=log(mean(i,1)-std(i)./sqrt(size(i,1)));
p1=plot(time, y2,'LineStyle','-', 'LineWidth',1.5,'Color',color_b)
patch([time fliplr(time)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time fliplr(time)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
hold on
box ('off')
%

for m=1:size(BZ.idrat,1)
    r=BZ.idrat(m);
new{m,1}=(BZ.power_thal{r,:});
end
z=vertcat(BZ.power_thal{:}); color_b= [0.5 0 0];
% z=BZ.power_ctx;color_b=[0.5 0.5 0.5];

y2=log(mean(z,1));
y1=log(mean(z,1)+std(z)./sqrt(size(z,1)));
y3=log(mean(z,1)-std(z)./sqrt(size(z,1)));
p2=plot(time, y2,'LineStyle','--', 'LineWidth',1.5,'Color',color_b)
patch([time fliplr(time)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time fliplr(time)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')

xlim ([0 80])
xticks([0:10:100])
legend([p1 p2],{'SNR','BZ'},'Box','off')
ylabel ('Subcortical Log Power')
xlabel ('Frequency(Hz)')
box ('off')
box ('off')
fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 10, 10];
fig.Color='w';



iz=[i;z];
fig=figure();
y2=log(mean(iz,1));
y1=log(mean(iz,1)+std(iz)./sqrt(size(iz,1)));
y3=log(mean(iz,1)-std(iz)./sqrt(size(iz,1)));
p1=plot(time, y2,'LineStyle','-', 'LineWidth',1.5,'Color',color_b)
% legend([p1],{'M1'},'Box','off')
patch([time fliplr(time)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time fliplr(time)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
hold on
box ('off')
xlim ([0 80])
xticks([0:10:100])
ylabel ('Cortical Log Power')
xlabel ('Frequency(Hz)')
box ('off')
box ('off')
fig.Units   = 'centimeters';
fig.OuterPosition= [10, 10, 10, 10];
fig.Color='w';
set(gca,'FontSize',12)


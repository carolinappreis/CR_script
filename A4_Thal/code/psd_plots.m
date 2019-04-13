clear all
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal')
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal')
load('BZ.mat');  
load('SNR.mat'); 
color_s=[0 0 0];
color_ctxb=[0.5 0.5 0.5];
time=1:501;

% i=SNR.power_ctx;
for m=1:size(SNR.idrat,1)
    r=SNR.idrat(m);
new{m,1}=(SNR.power_thal{r,:});
end
i=vertcat(new{:}); 
color_b= [0 0 0.5];
y2=log(mean(i,1));
y1=log(mean(i,1)+std(i)./sqrt(size(i,1)));
y3=log(mean(i,1)-std(i)./sqrt(size(i,1)));
p1=plot(time, y2,'LineStyle','-', 'LineWidth',1.5,'Color',color_b)
patch([time fliplr(time)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time fliplr(time)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
hold on
box ('off')
%
% z=BZ.power_ctx;
for m=1:size(BZ.idrat,1)
    r=BZ.idrat(m);
new{m,1}=(BZ.power_thal{r,:});
end
z=vertcat(BZ.power_thal{:}); color_b= [0.5 0 0.5];

y2=log(mean(z,1));
y1=log(mean(z,1)+std(z)./sqrt(size(z,1)));
y3=log(mean(z,1)-std(z)./sqrt(size(z,1)));
p2=plot(time, y2,'LineStyle','-', 'LineWidth',1.5,'Color',color_b)
patch([time fliplr(time)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time fliplr(time)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')

xlim ([0 80])
xticks([0:10:100])
legend([p1 p2],{'SNR','BZ'})
ylabel ('Subcortical Log Power')
xlabel ('Frequency(Hz)')
box ('off')



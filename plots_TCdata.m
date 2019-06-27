clear all
close all

%---------------------------------SNR------------------------------
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal')
%  cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal')
load('SNR_opt.mat');
time=1:1001;
time2=1:401;
for q=1:size (SNR.across_b_long,1)
    b_l_rat(q,:)=mean(SNR.across_b_long{q,1},1);
    b_s_rat(q,:)=mean(SNR.across_b_short{q,1},1);
    t_l_rat(q,:)=mean(SNR.across_t_long{q,1},1);
    t_s_rat(q,:)=mean(SNR.across_t_short{q,1},1);
    b_surr_rat(q,:)=mean(SNR.across_b_surr{q,1},1);
    t_surr_rat(q,:)=mean(SNR.across_t_surr{q,1},1);
end


%%% ------ SNR cortical burst
fig=figure()
color_b= [0 0 0.5]; 
y2=100.*(mean(SNR.pchange_ctxl)); 
y1=100.*(mean(SNR.pchange_ctxl)+std(SNR.pchange_ctxl)./sqrt(size(SNR.pchange_ctxl,1))); 
y3=100.*(mean(SNR.pchange_ctxl)-std(SNR.pchange_ctxl)./sqrt(size(SNR.pchange_ctxl,1)));
% y1= 100.*(mean(SNR.pchange_ctxl)+std(SNR.pchange_ctxl));
% y3= 100.*(mean(SNR.pchange_ctxl)-std(SNR.pchange_ctxl));
p2=plot(time, y2,'LineStyle','-', 'LineWidth',1.5,'Color',color_b)
patch([time fliplr(time)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time fliplr(time)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
xline(500,'--',{'burst onset'},'LabelOrientation','horizontal','Color',[0.5 0.5 0.5],'LineWidth',2)
hold on

color_b= [0.0118    0.7922    0.9882]; 
y2=100.*(mean(SNR.pchange_ctxs)); 
y1=100.*(mean(SNR.pchange_ctxs)+std(SNR.pchange_ctxs)./sqrt(size(SNR.pchange_ctxs,1))); 
y3=100.*(mean(SNR.pchange_ctxs)-std(SNR.pchange_ctxs)./sqrt(size(SNR.pchange_ctxs,1)));
% y1= 100.*(mean(SNR.pchange_ctxs)+std(SNR.pchange_ctxs));
% y3= 100.*(mean(SNR.pchange_ctxs)-std(SNR.pchange_ctxs));
p1=plot(time, y2,'LineStyle','-', 'LineWidth',1.5,'Color',color_b)
patch([time fliplr(time)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time fliplr(time)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
xline(500,'--',{'burst onset'},'LabelOrientation','horizontal','Color',[0.5 0.5 0.5],'LineWidth',2)

box('off')
xlim ([0 1000])
ylim ([-150 150])
xticks([0:250:1000])
xticklabels ({'-500','-250','0','250','500'})
fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 10, 10];
fig.Color='w';
set(gca,'FontSize',12)
ylabel('Beta amplitude change (%)');
xlabel('Time(msec)');
%legend([p1 p2],{'short burst','long burst'},'Box','off','Location','southeast')

%%% ------ SNR beta activity change

fig=figure()
st=NaN(1,1001);
clear A; A=SNR.pchange_l; 
clear B; B=SNR.pchange_surr;
hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
beg=find(st(1,:)<0.05 & st2(1,:)~=0);
if ~isempty(beg)
    sig_rise_all=[beg(1) beg(end)]; 
end
clear st st2

color_b= [0 0 0.5]; 
y2=100.*(mean(SNR.pchange_l)); 
y1=100.*(mean(SNR.pchange_l)+std(SNR.pchange_l)./sqrt(size(SNR.pchange_l,1))); 
y3=100.*(mean(SNR.pchange_l)-std(SNR.pchange_l)./sqrt(size(SNR.pchange_l,1)));
% y1= 100.*(mean(SNR.pchange_l)+std(SNR.pchange_l));
% y3= 100.*(mean(SNR.pchange_l)-std(SNR.pchange_l));
p2=plot(time, y2,'LineStyle','-', 'LineWidth',1.5,'Color',color_b)
patch([time fliplr(time)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time fliplr(time)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
if ~isempty (beg)
patch([sig_rise_all(1) sig_rise_all(2) sig_rise_all(2) sig_rise_all(1)],[15 15 15.5 15.5],color_b,'EdgeColor','none')
clear beg
end
xline(500,'--',{'burst onset'},'LabelOrientation','horizontal','Color',[0.5 0.5 0.5],'LineWidth',2)

hold on

st=NaN(1,1001);
clear A; A=SNR.pchange_s; 
clear B; B=SNR.pchange_surr;
hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
beg=find(st(1,:)<0.05 & st2(1,:)~=0);
if ~isempty(beg)
    sig_rise_all=[beg(1) beg(end)]; 
end
clear st st2

color_b= [0.0118    0.7922    0.9882]; 
y2=100.*(mean(SNR.pchange_s)); 
y1=100.*(mean(SNR.pchange_s)+std(SNR.pchange_s)./sqrt(size(SNR.pchange_s,1))); 
y3=100.*(mean(SNR.pchange_s)-std(SNR.pchange_s)./sqrt(size(SNR.pchange_s,1)));
% y1= 100.*(mean(SNR.pchange_s)+std(SNR.pchange_s));
% y3= 100.*(mean(SNR.pchange_s)-std(SNR.pchange_s));
p1=plot(time, y2,'LineStyle','-', 'LineWidth',1.5,'Color',color_b)
patch([time fliplr(time)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time fliplr(time)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
if ~isempty (beg)
patch([sig_rise_all(1) sig_rise_all(2) sig_rise_all(2) sig_rise_all(1)],[16 16 16.5 16.5],color_b,'EdgeColor','none')
clear beg
end

box('off')
xlim ([0 1000])
ylim ([-20 20])
xticks([0:250:1000])
xticklabels ({'-500','-250','0','250','500'})
fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 10, 10];
fig.Color='w';
set(gca,'FontSize',12)
ylabel('Beta amplitude change (%)');
xlabel('Time(msec)');
%legend([p1 p2],{'short burst','long burst'},'Box','off','Location','southeast')

%%%------SNR phase coupling across bursts
fig=figure()
color_b= [0 0 0.5]; 

st=NaN(1,401);
clear A; A=t_l_rat; 
clear B; B=t_surr_rat; 
hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
beg=find(st(1,:)<0.05 & st2(1,:)~=0);
if ~isempty(beg)
    sig_rise_all=[beg(1) beg(end)]; 
end
clear st st2

y2=zscore(mean(t_l_rat)); 
y1=zscore(mean(t_l_rat)+std(t_l_rat)./sqrt(size(t_l_rat,1))); 
y3=zscore(mean(t_l_rat)-std(t_l_rat)./sqrt(size(t_l_rat,1)));
p2=plot(time2, y2,'LineStyle','-', 'LineWidth',1.5,'Color',color_b)
patch([time2 fliplr(time2)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time2 fliplr(time2)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
hold on
if ~isempty(beg)
patch([sig_rise_all(1) sig_rise_all(2) sig_rise_all(2) sig_rise_all(1)],[2 2 2.05 2.05],color_b,'EdgeColor','none')
clear beg
end

hold on

st=NaN(1,401);
clear A; A=t_s_rat; 
clear B; B=t_surr_rat; 
hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
beg=find(st(1,:)<0.05 & st2(1,:)~=0);
if ~isempty(beg)
    sig_rise_all=[beg(1) beg(end)]; 
end
clear st st2

color_b= [0.0118    0.7922    0.9882]; 
y2=zscore(mean(t_s_rat)); 
y1=zscore(mean(t_s_rat)+std(t_s_rat)./sqrt(size(t_s_rat,1))); 
y3=zscore(mean(t_s_rat)-std(t_s_rat)./sqrt(size(t_s_rat,1)));
p1=plot(time2, y2,'LineStyle','-', 'LineWidth',1.5,'Color',color_b)
patch([time2 fliplr(time2)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time2 fliplr(time2)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
hold on
if ~isempty(beg)
patch([sig_rise_all(1) sig_rise_all(2) sig_rise_all(2) sig_rise_all(1)],[2.1 2.1 2.15 2.15],color_b,'EdgeColor','none')
clear beg
end
xline(200,'--',{'burst onset'},'LabelOrientation','horizontal','Color',[0.5 0.5 0.5],'LineWidth',2)

box('off')
xlim ([0 400])
ylim ([-3 3])
xticks([0:100:400])
xticklabels ({'-200','-100','0','100','200'})
fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 10, 10];
fig.Color='w';
set(gca,'FontSize',12)
ylabel('PSI across {\beta} burts (zscore)');
xlabel('Time(msec)');
%legend([p1 p2],{'short burst','long burst'},'Box','off','Location','southeast')

%%%------SNR phase coupling within bursts
fig=figure()
color_b= [0 0 0.5]; 

st=NaN(1,401);
clear A; A=b_l_rat; 
clear B; B=b_surr_rat; 
hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
beg=find(st(1,:)<0.05 & st2(1,:)~=0);
if ~isempty(beg)
    sig_rise_all=[beg(1) beg(end)]; 
end
clear st st2

y2=zscore(mean(b_l_rat)); 
y1=zscore(mean(b_l_rat)+std(b_l_rat)./sqrt(size(b_l_rat,1))); 
y3=zscore(mean(b_l_rat)-std(b_l_rat)./sqrt(size(b_l_rat,1)));
p2=plot(time2, y2,'LineStyle','-', 'LineWidth',1.5,'Color',color_b)
patch([time2 fliplr(time2)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time2 fliplr(time2)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
hold on
if ~isempty(beg)
patch([sig_rise_all(1) sig_rise_all(2) sig_rise_all(2) sig_rise_all(1)],[2 2 2.05 2.05],color_b,'EdgeColor','none')
clear beg
end

hold on

st=NaN(1,401);
clear A; A=b_s_rat; 
clear B; B=b_surr_rat; 
hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
beg=find(st(1,:)<0.05 & st2(1,:)~=0);
if ~isempty(beg)
    sig_rise_all=[beg(1) beg(end)]; 
end
clear st st2

color_b= [0.0118    0.7922    0.9882]; 
y2=zscore(mean(b_s_rat)); 
y1=zscore(mean(b_s_rat)+std(b_s_rat)./sqrt(size(b_s_rat,1))); 
y3=zscore(mean(b_s_rat)-std(b_s_rat)./sqrt(size(b_s_rat,1)));
p1=plot(time2, y2,'LineStyle','-', 'LineWidth',1.5,'Color',color_b)
patch([time2 fliplr(time2)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time2 fliplr(time2)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
hold on
if ~isempty(beg)
patch([sig_rise_all(1) sig_rise_all(2) sig_rise_all(2) sig_rise_all(1)],[2.1 2.1 2.15 2.15],color_b,'EdgeColor','none')
clear beg
end
xline(200,'--',{'burst onset'},'LabelOrientation','horizontal','Color',[0.5 0.5 0.5],'LineWidth',2)

box('off')
xlim ([0 400])
ylim ([-3 3])
xticks([0:100:400])
xticklabels ({'-200','-100','0','100','200'})
fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 10, 10];
fig.Color='w';
set(gca,'FontSize',12)
ylabel('PSI within {\beta} burts (zscore)');
xlabel('Time(msec)');
%legend([p1 p2],{'short burst','long burst'},'Box','off','Location','southeast')





%% ---------------------------------BZ------------------------------
clear all
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal')
% cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal')
load('BZ_opt.mat');
time=1:1001;
time2=1:401;
for q=1:size (BZ.across_b_long,1)
    b_l_rat(q,:)=mean(BZ.across_b_long{q,1},1);
    b_s_rat(q,:)=mean(BZ.across_b_short{q,1},1);
    t_l_rat(q,:)=mean(BZ.across_t_long{q,1},1);
    t_s_rat(q,:)=mean(BZ.across_t_short{q,1},1);
    b_surr_rat(q,:)=mean(BZ.across_b_surr{q,1},1);
    t_surr_rat(q,:)=mean(BZ.across_t_surr{q,1},1);
end

%%% ------ BZ cortical burst

fig=figure()
color_b= [0.5 0 0]; 
y2=100.*(mean(BZ.pchange_ctxl)); 
y1=100.*mean(BZ.pchange_ctxl)+std(BZ.pchange_ctxl)./sqrt(size(BZ.pchange_ctxl,1)); 
y3=100.*mean(BZ.pchange_ctxl)-std(BZ.pchange_ctxl)./sqrt(size(BZ.pchange_ctxl,1));
% y1= 100.*(mean(BZ.pchange_ctxl)+std(BZ.pchange_ctxl));
% y3= 100.*(mean(BZ.pchange_ctxl)-std(BZ.pchange_ctxl));
p2=plot(time, y2,'LineStyle','-', 'LineWidth',1.5,'Color',color_b)
patch([time fliplr(time)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time fliplr(time)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
xline(500,'--',{'burst onset'},'LabelOrientation','horizontal','Color',[0.5 0.5 0.5],'LineWidth',2)
hold on

color_b= [0.7098    0.5725    0.0157]; 
y2=100.*(mean(BZ.pchange_ctxs)); 
y1=100.*(mean(BZ.pchange_ctxs)+std(BZ.pchange_ctxs)./sqrt(size(BZ.pchange_ctxs,1))); 
y3=100.*(mean(BZ.pchange_ctxs)-std(BZ.pchange_ctxs)./sqrt(size(BZ.pchange_ctxs,1)));
% y1= 100.*(mean(BZ.pchange_ctxs)+std(BZ.pchange_ctxs));
% y3= 100.*(mean(BZ.pchange_ctxs)-std(BZ.pchange_ctxs));
p1=plot(time, y2,'LineStyle','-', 'LineWidth',1.5,'Color',color_b)
patch([time fliplr(time)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time fliplr(time)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
xline(500,'--',{'burst onset'},'LabelOrientation','horizontal','Color',[0.5 0.5 0.5],'LineWidth',2)

box('off')
xlim ([0 1000])
ylim ([-150 150])
xticks([0:250:1000])
xticklabels ({'-500','-250','0','250','500'})
fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 10, 10];
fig.Color='w';
set(gca,'FontSize',12)
ylabel('Beta amplitude change (%)');
xlabel('Time(msec)');
%legend([p1 p2],{'short burst','long burst'},'Box','off','Location','southeast')

%%% ------ BZ beta activity

fig=figure()
st=NaN(1,1001);
clear A; A=BZ.pchange_l; 
clear B; B=BZ.pchange_surr;
hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
beg=find(st(1,:)<0.05 & st2(1,:)~=0);
if ~isempty(beg)
    sig_rise_all=[beg(1) beg(end)]; 
end
clear st st2

color_b= [0.5 0 0]; 
y2=100.*(mean(BZ.pchange_l)); 
y1=100.*(mean(BZ.pchange_l)+std(BZ.pchange_l)./sqrt(size(BZ.pchange_l,1))); 
y3=100.*(mean(BZ.pchange_l)-std(BZ.pchange_l)./sqrt(size(BZ.pchange_l,1)));
% y1= 100.*(mean(BZ.pchange_l)+std(BZ.pchange_l));
% y3= 100.*(mean(BZ.pchange_l)-std(BZ.pchange_l));
p2=plot(time, y2,'LineStyle','-', 'LineWidth',1.5,'Color',color_b)
patch([time fliplr(time)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time fliplr(time)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
if ~isempty (beg)
patch([sig_rise_all(1) sig_rise_all(2) sig_rise_all(2) sig_rise_all(1)],[15 15 15+0.5 15+0.5],color_b,'EdgeColor','none')
clear beg
end
xline(500,'--',{'burst onset'},'LabelOrientation','horizontal','Color',[0.5 0.5 0.5],'LineWidth',2)
hold on

st=NaN(1,1001);
clear A; A=BZ.pchange_s; 
clear B; B=BZ.pchange_surr;
hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
beg=find(st(1,:)<0.05 & st2(1,:)~=0);
if ~isempty(beg)
    sig_rise_all=[beg(1) beg(end)]; 
end
clear st st2

color_b= [0.7098    0.5725    0.0157]; 
y2=100.*(mean(BZ.pchange_s)); 
y1=100.*(mean(BZ.pchange_s)+std(BZ.pchange_s)./sqrt(size(BZ.pchange_s,1))); 
y3=100.*(mean(BZ.pchange_s)-std(BZ.pchange_s)./sqrt(size(BZ.pchange_s,1)));
% y1= 100.*(mean(BZ.pchange_s)+std(BZ.pchange_s));
% y3= 100.*(mean(BZ.pchange_s)-std(BZ.pchange_s));
p1=plot(time, y2,'LineStyle','-', 'LineWidth',1.5,'Color',color_b)
patch([time fliplr(time)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time fliplr(time)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
if ~isempty (beg)
patch([sig_rise_all(1) sig_rise_all(2) sig_rise_all(2) sig_rise_all(1)],[16 16 16.5 16.5],color_b,'EdgeColor','none')
clear beg
end
xline(500,'--',{'burst onset'},'LabelOrientation','horizontal','Color',[0.5 0.5 0.5],'LineWidth',2)

box('off')
xlim ([0 1000])
ylim ([-20 20])
xticks([0:250:1000])
xticklabels ({'-500','-250','0','250','500'})
fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 10, 10];
fig.Color='w';
set(gca,'FontSize',12)
ylabel('Beta amplitude change (%)');
xlabel('Time(msec)');
%legend([p1 p2],{'short burst','long burst'},'Box','off','Location','southeast')

%%%------BZ phase coupling across bursts
fig=figure()
color_b= [0.5 0 0]; 

st=NaN(1,401);
clear A; A=t_l_rat; 
clear B; B=t_surr_rat; 
hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
beg=find(st(1,:)<0.05 & st2(1,:)~=0);
if ~isempty(beg)
    sig_rise_all=[beg(1) beg(end)]; 
end
clear st st2

y2=zscore(mean(t_l_rat)); 
y1=zscore(mean(t_l_rat)+std(t_l_rat)./sqrt(size(t_l_rat,1))); 
y3=zscore(mean(t_l_rat)-std(t_l_rat)./sqrt(size(t_l_rat,1)));
p2=plot(time2, y2,'LineStyle','-', 'LineWidth',1.5,'Color',color_b)
patch([time2 fliplr(time2)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time2 fliplr(time2)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
hold on
if ~isempty(beg)
patch([sig_rise_all(1) sig_rise_all(2) sig_rise_all(2) sig_rise_all(1)],[2 2 2.05 2.05],color_b,'EdgeColor','none')
clear beg
end

hold on

st=NaN(1,401);
clear A; A=t_s_rat; 
clear B; B=t_surr_rat; 
hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
beg=find(st(1,:)<0.05 & st2(1,:)~=0);
if ~isempty(beg)
    sig_rise_all=[beg(1) beg(end)]; 
end
clear st st2

color_b= [0.7098    0.5725    0.0157]; 
y2=zscore(mean(t_s_rat)); 
y1=zscore(mean(t_s_rat)+std(t_s_rat)./sqrt(size(t_s_rat,1))); 
y3=zscore(mean(t_s_rat)-std(t_s_rat)./sqrt(size(t_s_rat,1)));
p1=plot(time2, y2,'LineStyle','-', 'LineWidth',1.5,'Color',color_b)
patch([time2 fliplr(time2)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time2 fliplr(time2)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
hold on
if ~isempty(beg)
patch([sig_rise_all(1) sig_rise_all(2) sig_rise_all(2) sig_rise_all(1)],[2.1 2.1 2.15 2.15],color_b,'EdgeColor','none')
clear beg
end
xline(200,'--',{'burst onset'},'LabelOrientation','horizontal','Color',[0.5 0.5 0.5],'LineWidth',2)

box('off')
xlim ([0 400])
ylim ([-3 3])
xticks([0:100:400])
xticklabels ({'-200','-100','0','100','200'})
fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 10, 10];
fig.Color='w';
set(gca,'FontSize',12)
ylabel('PSI across {\beta} burts (zscore)');
xlabel('Time(msec)');
%legend([p1 p2],{'short burst','long burst'},'Box','off','Location','southeast')

%%%------BZ phase coupling within bursts
fig=figure()
color_b= [0.5 0 0]; 

st=NaN(1,401);
clear A; A=b_l_rat; 
clear B; B=b_surr_rat; 
hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
beg=find(st(1,:)<0.05 & st2(1,:)~=0);
if ~isempty(beg)
    sig_rise_all=[beg(1) beg(end)]; 
end
clear st st2

y2=zscore(mean(b_l_rat)); 
y1=zscore(mean(b_l_rat)+std(b_l_rat)./sqrt(size(b_l_rat,1))); 
y3=zscore(mean(b_l_rat)-std(b_l_rat)./sqrt(size(b_l_rat,1)));
p2=plot(time2, y2,'LineStyle','-', 'LineWidth',1.5,'Color',color_b)
patch([time2 fliplr(time2)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time2 fliplr(time2)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
hold on
if ~isempty(beg)
patch([sig_rise_all(1) sig_rise_all(2) sig_rise_all(2) sig_rise_all(1)],[2 2 2.05 2.05],color_b,'EdgeColor','none')
clear beg
end

hold on

st=NaN(1,401);
clear A; A=b_s_rat; 
clear B; B=b_surr_rat; 
hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
beg=find(st(1,:)<0.05 & st2(1,:)~=0);
if ~isempty(beg)
    sig_rise_all=[beg(1) beg(end)]; 
end
clear st st2

color_b= [0.7098    0.5725    0.0157];
y2=zscore(mean(b_s_rat)); 
y1=zscore(mean(b_s_rat)+std(b_s_rat)./sqrt(size(b_s_rat,1))); 
y3=zscore(mean(b_s_rat)-std(b_s_rat)./sqrt(size(b_s_rat,1)));
p1=plot(time2, y2,'LineStyle','-', 'LineWidth',1.5,'Color',color_b)
patch([time2 fliplr(time2)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time2 fliplr(time2)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
hold on
if ~isempty(beg)
patch([sig_rise_all(1) sig_rise_all(2) sig_rise_all(2) sig_rise_all(1)],[2.1 2.1 2.15 2.15],color_b,'EdgeColor','none')
clear beg
end
xline(200,'--',{'burst onset'},'LabelOrientation','horizontal','Color',[0.5 0.5 0.5],'LineWidth',2)

box('off')
xlim ([0 400])
ylim ([-3 3])
xticks([0:100:400])
xticklabels ({'-200','-100','0','100','200'})
fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 10, 10];
fig.Color='w';
set(gca,'FontSize',12)
ylabel('PSI within {\beta} burts (zscore)');
xlabel('Time(msec)');
%legend([p1 p2],{'short burst','long burst'},'Box','off','Location','southeast')


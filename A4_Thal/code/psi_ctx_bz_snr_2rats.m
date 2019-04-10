clear all
% cd ('C:\Users\creis\Documents\GitHub\CRcode\codes_thal\A4_Thal')
cd('/Users/Carolina/Documents/GitHub/CRcode/codes_thal/A4_Thal')
load 'bz_rat_level.mat'; color_b=[0.5 0 0.5];


color_s=[0 0 0];
color_ctxb=[0.5 0.5 0.5];
time= [1:2*200+1];

reg_m=mean(clust_b_m(end-1:end,:));
reg_sd=std(clust_b_m(end-1:end,:))./sqrt(2);
reg_sm=mean(clust_s_m(end-1:end,:));
reg_ssd=std(clust_s_m(end-1:end,:))./sqrt(2);


cd ('/Users/Carolina/Documents/GitHub/CRcode/codes_thal/A4_Thal/code')
st=NaN(1,401);
clear A; A=clust_b_m(end-1:end,:)
clear B; B=clust_s_m(end-1:end,:)
hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
beg=find(st(1,:)<0.01 & st2(1,:)~=0);
if ~isempty(beg)
    sig_rise_all=[beg(1) beg(end)];
    
end

y2=reg_m; y1=y2+reg_sd; y3=y2-reg_sd;
y5=reg_sm; y4=y5+reg_ssd; y6=y5-reg_ssd;
p1=plot(time, y2,'LineStyle','-.', 'LineWidth',1.5,'Color',color_b)
patch([time fliplr(time)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time fliplr(time)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
hold on
p2=plot(time, y5,'LineStyle','-.','LineWidth',1.5,'Color',color_s);
patch([time fliplr(time)], [y4 fliplr(y5)],[color_s],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time fliplr(time)], [y5 fliplr(y6)],[color_s],'FaceAlpha',[0.2],'EdgeColor','none')
xlim ([0 400])
ylim ([0 0.6])
xticks([0:100:400])
xticklabels ({'-200','-100','0','100','200'})

patch([sig_rise_all(1) sig_rise_all(2) sig_rise_all(2) sig_rise_all(1)],[-0.1 -0.1 0.6 0.6],color_b,'FaceAlpha',[0.1],'EdgeColor','none')

hold on

cd('/Users/Carolina/Documents/GitHub/CRcode/codes_thal/A4_Thal')
load 'snr_rat_level.mat'; color_b=[0 0 0.5];

reg_m=mean(clust_b_m(end-1:end,:));
reg_sd=std(clust_b_m(end-1:end,:))./sqrt(2);
reg_sm=mean(clust_s_m(end-1:end,:));
reg_ssd=std(clust_s_m(end-1:end,:))./sqrt(2);


cd ('/Users/Carolina/Documents/GitHub/CRcode/codes_thal/A4_Thal/code')
st=NaN(1,401);
clear A; A=clust_b_m(end-1:end,:)
clear B; B=clust_s_m(end-1:end,:)
hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
beg=find(st(1,:)<0.01 & st2(1,:)~=0);
if ~isempty(beg)
    sig_rise_all=[beg(1) beg(end)];
    
end


y2=reg_m; y1=y2+reg_sd; y3=y2-reg_sd;
y5=reg_sm; y4=y5+reg_ssd; y6=y5-reg_ssd;
p3=plot(time, y2, 'LineWidth',1.5,'Color',color_b)
patch([time fliplr(time)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time fliplr(time)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
hold on
p4=plot(time, y5,'LineWidth',1.5,'Color',color_s);
patch([time fliplr(time)], [y4 fliplr(y5)],[color_s],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time fliplr(time)], [y5 fliplr(y6)],[color_s],'FaceAlpha',[0.2],'EdgeColor','none')
xlim ([0 400])
ylim ([0 0.6])
xticks([0:100:400])
xticklabels ({'-200','-100','0','100','200'})

patch([sig_rise_all(1) sig_rise_all(2) sig_rise_all(2) sig_rise_all(1)],[-0.1 -0.1 0.6 0.6],color_b,'FaceAlpha',[0.1],'EdgeColor','none')

legend([p1 p2 p3 p4],{'BZ-CTX in burst','BZ-CTX Baseline','SNr-CTX in burst','SNr-CTX Baseline'})
ylabel ('Phase Synchrony Index')
xlabel ('Time (msec)')
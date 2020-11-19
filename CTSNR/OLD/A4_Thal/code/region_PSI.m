clear all
cd ('C:\Users\creis\Documents\GitHub\CRcode\codes_thal\A4_Thal')
% cd('/Users/Carolina/Documents/GitHub/CRcode/codes_thal/A4_Thal')
load 'bz_allrats.mat'; color_b=[0.5 0 0.5];id_rat=([1 7 9 11])'; rats=size(id_rat,1);
% id_rat=([randi([1,5],1,1) randi([6,7],1,1) randi([8,10],1,1) 11]);
color_s=[0 0 0];
color_ctxb=[0.5 0.5 0.5];
time= [1:2*200+1];
 
%-----

for i =1:length(id_rat)
clust_b{id_rat(i),:};
clust_b_m(i,:)=mean(ans);
clust_b_sd(i,:)=std(ans)./sqrt(size(ans,1));
clust_s{id_rat(i),:};
clust_s_m(i,:)=mean(ans);
clust_s_sd(i,:)=std(ans)./sqrt(size(ans,1));
power_coh{id_rat(i),:};
power_coh_m(i,:)=mean(ans);
power_coh_sd(i,:)=std(ans)./sqrt(size(ans,1));
power_ctx1(i,:)=power_ctx{id_rat(i),1}(1,:);
end

ctx_b=ctx_b1(id_rat,:);
ctx_b_m=mean(ctx_b,1);
ctx_b_sd=std(ctx_b)./sqrt(size(ctx_b,1));

power_ctx_m=mean(power_ctx1);
power_ctx_sd=std(power_ctx1)./sqrt(size(power_ctx1,1));

%-----

reg_m=mean(clust_b_m);
reg_sd=std(clust_b_m)./sqrt(size(clust_b_m,1));
reg_sm=mean(clust_s_m);
reg_ssd=std(clust_s_m)./sqrt(size(clust_s_m,1));


% cd ('/Users/Carolina/Documents/GitHub/CRcode/codes_thal/A4_Thal/code')
cd('C:\Users\creis\Documents\GitHub\CRcode\codes_thal\A4_Thal\code')
st=NaN(1,401);
clear A; A=clust_b_m; %b1{f,1};
clear B; B=clust_s_m; %s1{f,1}(1:size(A,1),:);
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

patch([sig_rise_all(1) sig_rise_all(2) sig_rise_all(2) sig_rise_all(1)],[0 0 1 1],color_b,'FaceAlpha',[0.1],'EdgeColor','none')

hold on

clearvars -except p1 p2
cd ('C:\Users\creis\Documents\GitHub\CRcode\codes_thal\A4_Thal')
% cd('/Users/Carolina/Documents/GitHub/CRcode/codes_thal/A4_Thal')
load 'snr_allrats.mat'; color_b=[0 0 0.5];id_rat=([1 3 11 12])'; rats=size(id_rat,1);
% id_rat=([randi([1,2],1,1) randi([3,5],1,1) randi([6,11],1,1) randi([12,14],1,1)]);
color_s=[0 0 0];
color_ctxb=[0.5 0.5 0.5];
time= [1:2*200+1];

%-----

for i =1:length(id_rat)
clust_b{id_rat(i),:};
clust_b_m(i,:)=mean(ans);
clust_b_sd(i,:)=std(ans)./sqrt(size(ans,1));
clust_s{id_rat(i),:};
clust_s_m(i,:)=mean(ans);
clust_s_sd(i,:)=std(ans)./sqrt(size(ans,1));
power_coh{id_rat(i),:};
power_coh_m(i,:)=mean(ans);
power_coh_sd(i,:)=std(ans)./sqrt(size(ans,1));
power_ctx1(i,:)=power_ctx{id_rat(i),1}(1,:);
end

ctx_b=ctx_b1(id_rat,:);
ctx_b_m=mean(ctx_b,1);
ctx_b_sd=std(ctx_b)./sqrt(size(ctx_b,1));

power_ctx_m=mean(power_ctx1);
power_ctx_sd=std(power_ctx1)./sqrt(size(power_ctx1,1));

%-----

reg_m=mean(clust_b_m);
reg_sd=std(clust_b_m)./sqrt(size(clust_b_m,1));
reg_sm=mean(clust_s_m);
reg_ssd=std(clust_s_m)./sqrt(size(clust_s_m,1));


% cd ('/Users/Carolina/Documents/GitHub/CRcode/codes_thal/A4_Thal/code')
cd('C:\Users\creis\Documents\GitHub\CRcode\codes_thal\A4_Thal\code')
st=NaN(1,401);
clear A; A=clust_b_m; %b1{f,1};
clear B; B=clust_s_m; %s1{f,1}(1:size(A,1),:);
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

patch([sig_rise_all(1) sig_rise_all(2) sig_rise_all(2) sig_rise_all(1)],[0 0 1 1],color_b,'FaceAlpha',[0.1],'EdgeColor','none')

legend([p1 p2 p3 p4],{'BZ-CTX in burst','BZ-CTX Baseline','SNr-CTX in burst','SNr-CTX Baseline'})
ylabel ('Phase Synchrony Index')
xlabel ('Time (msec)')
ylim([0 1])
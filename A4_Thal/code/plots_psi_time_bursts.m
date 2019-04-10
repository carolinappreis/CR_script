clear all
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal')
% load('BZ.mat'); color_b= [0.5 0 0.5]; color_b1=color_b+[0.3 0 0.3];
load('SNR.mat'); color_b= [0 0 0.5]; color_b1=color_b+[0 0 0.3];
time=1:401;

for q=1:size (SNR.across_b_long,1)
    b_l_rat(q,:)=mean(SNR.across_b_long{q,1},1);
    b_s_rat(q,:)=mean(SNR.across_b_short{q,1},1);
    t_l_rat(q,:)=mean(SNR.across_t_long{q,1},1);
    t_s_rat(q,:)=mean(SNR.across_t_short{q,1},1);
    b_surr_rat(q,:)=mean(SNR.across_b_surr{q,1},1);
    t_surr_rat(q,:)=mean(SNR.across_t_surr{q,1},1);
end
% 
% for q=1:size (SNR.across_b_long,1)
%     b_psl_rat(q,:)=mean(SNR.phashift_b_long{q,1},1);
%     b_pss_rat(q,:)=mean(SNR.phashift_b_short{q,1},1);
%     t_psl_rat(q,:)=mean(SNR.phashift_t_long{q,1},1);
%     t_pss_rat(q,:)=mean(SNR.phashift_t_short{q,1},1);
% end


cd('C:\Users\creis\Documents\GitHub\CR_script\A4_Thal\code')
st=NaN(1,401);
clear A; A=t_s_rat; %b1{f,1};
clear B; B=t_surr_rat; %s1{f,1}(1:size(A,1),:);
hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
beg=find(st(1,:)<0.01 & st2(1,:)~=0);
if ~isempty(beg)
    sig_rise_all=[beg(1) beg(end)]; 
end

y2=zscore(mean(t_s_rat)); 
y1=zscore(mean(t_s_rat)+std(t_s_rat)./sqrt(size(t_s_rat,1))); 
y3=zscore(mean(t_s_rat)-std(t_s_rat)./sqrt(size(t_s_rat,1)));
p1=plot(time, y2,'LineStyle','-', 'LineWidth',1.5,'Color',color_b1)
patch([time fliplr(time)], [y1 fliplr(y2)],[color_b1],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time fliplr(time)], [y2 fliplr(y3)],[color_b1],'FaceAlpha',[0.2],'EdgeColor','none')
hold on
xlim ([0 400])
ylim ([-3 3])
xticks([0:100:400])
xticklabels ({'-200','-100','0','100','200'})
patch([sig_rise_all(1) sig_rise_all(2) sig_rise_all(2) sig_rise_all(1)],[2.9 2.9 3 3],color_b1,'EdgeColor','none')

hold on

cd('C:\Users\creis\Documents\GitHub\CR_script\A4_Thal\code')
st=NaN(1,401);
clear A; A=t_l_rat; %b1{f,1};
clear B; B=t_surr_rat; %s1{f,1}(1:size(A,1),:);
hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
beg=find(st(1,:)<0.01 & st2(1,:)~=0);
if ~isempty(beg)
    sig_rise_all=[beg(1) beg(end)]; 
end

y2=zscore(mean(t_l_rat)); 
y1=zscore(mean(t_l_rat)+std(t_l_rat)./sqrt(size(t_l_rat,1))); 
y3=zscore(mean(t_l_rat)-std(t_l_rat)./sqrt(size(t_l_rat,1)));
p2=plot(time, y2,'LineStyle','-', 'LineWidth',1.5,'Color',color_b)
patch([time fliplr(time)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time fliplr(time)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
hold on
patch([sig_rise_all(1) sig_rise_all(2) sig_rise_all(2) sig_rise_all(1)],[2.5 2.5 2.6 2.6],color_b,'EdgeColor','none')

legend([p1 p2],{'short burst','long burst'})
ylabel ('PSI bursts (zscore)')
xlabel ('Time (msec)')


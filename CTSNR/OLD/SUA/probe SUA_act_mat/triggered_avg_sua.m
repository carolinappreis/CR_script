clear all
%   cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\probe SUA_act_mat')
 cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/probe SUA_act_mat')
load('BZ_sua_skrate.mat') ; 
subj= BZ.beta_rats;
% clear all
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\probe SUA_act_mat')
% load('BZ_sua_skrate.mat') ;
% subj= BZ.beta_rats;

Ecog=BZ.ctx(subj,:);
for i =1:size(subj,1)
    for ii=1:size(BZ.beta_idx{subj(i),1},2)
        hp=BZ.beta_idx{subj(i),1}(1,ii);
        data_all{i,1}(ii,:)=BZ.sua{subj(i),1}(hp,:);
        clear hp
    end
end

srn=1000;

% % % check that all units firing at beta freqs
% % data_a=cell2mat(data_all);
% % for ii=1:size(data_a,1)
% %     data=data_a(ii,:);
% %     
% %     spkrate_1=[];
% %     for i =1:srn:(length(data)-srn);
% %         spkrate_1=[spkrate_1 numel(find(data(i:i+srn)==1))];
% %     end
% %     srate(ii,1)=mean(spkrate_1);
% %     clear data
% % end

for j=1:size(data_all,1)
    
    data=data_all{j,1};
    ecog=Ecog(j,:);
    clearvars -except j ii BZ Ecog data ecog data_all srn rec_pa1 rec_npa1 rec_pa2 rec_npa2 region_pl region_spl region_npl region_nspl subj color_b
    
    for ii=1:size(data,1)
        
        [b,a]=butter(2,[15/(0.5*srn) 35/(0.5*srn)],'bandpass');
        Ecogfiltered=filtfilt(b,a,ecog);
        env=abs(hilbert(Ecogfiltered));
        onset1=bursts(env);
        onset1=horzcat(onset1{:});
        onset=bursts_aligned(env,Ecogfiltered);
        onset=horzcat(onset{:});
        
        data_g=smoothdata(data(ii,:),'gaussian',20);
        
        for jj=1:size(onset,2)
            if onset(jj)>200 && onset(jj)+200<length(data_g)
                output_pa(ii,jj,:)=data_g(onset(jj)-200:onset(jj)+200);
                output_npa(ii,jj,:)= data_g(onset1(jj)-200:onset1(jj)+200);
                output_count(ii,jj,:)= data(onset1(jj)-200:onset1(jj)+200);
            end
        end
        rec_pa1{j,1}=reshape(sum(output_pa,2),(size(sum(output_pa,2),1)),(size(sum(output_pa,2),3)));
        rec_npa1{j,1}=reshape(sum(output_npa,2),(size(sum(output_npa,2),1)),(size(sum(output_npa,2),3)));
    end
    rec_pa2(j,:)=mean(rec_pa1{j,1},1);
    rec_npa2(j,:)=mean(rec_npa1{j,1},1);
    clear data
end

rec_pa=cell2mat(rec_pa1);
rec_npa=cell2mat(rec_npa1);

clear j
for j=1:size(rec_pa,1)
    filt_pa2(j,:)=filtfilt(b,a,rec_pa(j,:));
    filt_npa2(j,:)=filtfilt(b,a,rec_npa(j,:));
    env_pa2(j,:)=abs(hilbert(filt_pa2(j,:)));
    env_npa2(j,:)=abs(hilbert(filt_npa2(j,:)));
    stat_pa2(j,1:2)=[sum(env_pa2(j,1:200)) sum(env_pa2(j,201:400))];
    stat_npa2(j,1:2)=[sum(env_npa2(j,1:200)) sum(env_npa2(j,201:400))];
end

stackedplot(rec_pa')
stackedplot(filt_pa2')

[p(1) h(1)]=ttest(stat_pa2(:,1),stat_pa2(:,2));
[p(2) h(2)]=ttest(stat_npa2(:,1),stat_npa2(:,2));
[p(3) h(3)]=ttest(stat_pa2(:,1),stat_npa2(:,1));
[p(4) h(4)]=ttest(stat_pa2(:,2),stat_npa2(:,2));

figure;
subplot(1,2,1)
plot(env_pa2')
hold on
plot(mean(env_pa2),'k','LineWidth',2)
ylim([0 0.8])
xlim([0 400])
box('off')
title('phase-aligned envelope')


subplot(1,2,2)
plot(env_npa2')
hold on
plot(mean(env_npa2),'k','LineWidth',2)
ylim([0 0.8])
xlim([0 400])
box('off')
title('non-phase-aligned envelope')


fig=figure;
subplot(1,2,1)
bar(1:2,[median(stat_pa2(:,1)) median(stat_pa2(:,2))],'FaceColor',[0.5 0.5 0.5],'FaceAlpha',.7,'EdgeColor','none','BarWidth',0.5);
hold on
plot(1,stat_pa2(:,1),'.')
plot(2,stat_pa2(:,2),'.')
box('off')
xlim([0.5 2.5])
title('phase-aligned envelope')
txt=(sprintf('p= %d',( round(h(1),4))));
text(1,58,txt)
set(gca,'FontSize',12)

subplot(1,2,2)

bar(1:2,[median(stat_npa2(:,1)) median(stat_npa2(:,2))],'FaceColor',[0.5 0.5 0.5],'FaceAlpha',.7,'EdgeColor','none','BarWidth',0.5);
hold on
plot(1,stat_npa2(:,1),'.')
plot(2,stat_npa2(:,2),'.')
box('off')
xlim([0.5 2.5])
title('non-phase-aligned envelope')
txt=(sprintf('p= %d',( round(h(2),4))));
text(1,58,txt)
set(gca,'FontSize',12)

fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 18, 15];
fig.Color='w';


region_pl=zscore(mean(rec_pa,1));
region_spl=zscore(std(rec_pa)./sqrt(size(rec_pa,1)));
region_npl=zscore(mean(rec_npa,1));
region_snpl=zscore(std(rec_npa)./sqrt(size(rec_npa,1)));


%%% just plot
time2=[1:401];
color_s=[0.5 0.5 0.5];
% color_b=[0 0 0.5]; %BZ
  color_b=[0.5 0 0]; %bz

fig=figure;
subplot(1,2,1)
y2=region_pl; y1=y2+region_spl; y3=y2-region_spl;
y5=region_npl; y4=y5+region_snpl; y6=y5-region_snpl;
p1=plot(time2, y2, 'LineWidth',1.5,'Color',color_b)
patch([time2 fliplr(time2)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.1],'EdgeColor','none')
patch([time2 fliplr(time2)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.1],'EdgeColor','none')
xline(200,'--',{'burst onset'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','Color',[0.5 0.5 0.5],'LineWidth',2)
xlim ([0 400])
ylim ([-5 5])
xticks([0:100:400])
xticklabels ({'-200','-100','0','100','200'})
legend([p1],{'phase-aligned'},'FontSize',12,'box','off','Location','northeast' )
box ('off')
xlabel ('Time (msec)')
ylabel('Firing-rate(z-score)')
set(gca,'FontSize',12)


subplot(1,2,2)
p2=plot(time2, y5, 'LineWidth',1.5,'Color',color_s)
patch([time2 fliplr(time2)], [y4 fliplr(y5)],[color_s],'FaceAlpha',[0.1],'EdgeColor','none')
patch([time2 fliplr(time2)], [y5 fliplr(y6)],[color_s],'FaceAlpha',[0.1],'EdgeColor','none')
xline(200,'--',{'burst onset'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','Color',[0.5 0.5 0.5],'LineWidth',2)
xlim ([0 400])
ylim ([-5 5])
xticks([0:100:400])
xticklabels ({'-200','-100','0','100','200'})
legend([p2],{'non-phase aligned'},'FontSize',12,'box','off','Location','northeast' )
box ('off')
xlabel ('Time (msec)')
ylabel('Firing-rate(z-score)')
set(gca,'FontSize',12)

fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 15, 10];
fig.Color='w';


[b,a]=butter(2,[20/(0.5*srn) 30/(0.5*srn)],'bandpass');

fig=figure;
subplot(1,2,1)
y2=filtfilt(b,a,region_pl); y1=y2+filtfilt(b,a,region_spl); y3=y2-filtfilt(b,a,region_spl);
y5=filtfilt(b,a,region_npl); y4=y5+filtfilt(b,a,region_snpl); y6=y5-filtfilt(b,a,region_spl);
p1=plot(time2, y2, 'LineWidth',1.5,'Color',color_b)
patch([time2 fliplr(time2)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.1],'EdgeColor','none')
patch([time2 fliplr(time2)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.1],'EdgeColor','none')
xline(200,'--',{'burst onset'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','Color',[0.5 0.5 0.5],'LineWidth',2)
xlim ([0 400])
ylim ([-5 5])
xticks([0:100:400])
xticklabels ({'-200','-100','0','100','200'})
legend([p1],{'phase-aligned'},'FontSize',12,'box','off','Location','northeast' )
box ('off')
xlabel ('Time (msec)')
ylabel('Firing-rate(z-score)')
set(gca,'FontSize',12)


subplot(1,2,2)
p2=plot(time2, y5, 'LineWidth',1.5,'Color',color_s)
patch([time2 fliplr(time2)], [y4 fliplr(y5)],[color_s],'FaceAlpha',[0.1],'EdgeColor','none')
patch([time2 fliplr(time2)], [y5 fliplr(y6)],[color_s],'FaceAlpha',[0.1],'EdgeColor','none')
xline(200,'--',{'burst onset'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','Color',[0.5 0.5 0.5],'LineWidth',2)
xlim ([0 400])
ylim ([-5 5])
xticks([0:100:400])
xticklabels ({'-200','-100','0','100','200'})
legend([p2],{'non-phase aligned'},'FontSize',12,'box','off','Location','northeast' )
box ('off')
xlabel ('Time (msec)')
ylabel('Firing-rate(z-score)')
set(gca,'FontSize',12)

fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 15, 10];
fig.Color='w';
% clear all
% cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/probe SUA_act_mat')
% % % load('BZ_rec_pa12.mat')
% b=rec_pa2;
% clear rec_pa2
% % % load('BZ_rec_pa12.mat')
% s=rec_pa2;
% clear rec_pa2
% for i=1:22
%     subplot(4,6,i)
%     plot(b(i,:))
%     hold on
%     plot(s(i,:))
% box('off')
% end



% for nn=3
% for i =1:size(s{nn,1},1)
% subplot(size(s{nn,1},1),1,i)
% plot(s{1,1}(i,:))
% xlim ([0 10000])
% end
% figure(2)
% for i =1:size(b{nn,1},1)
% subplot(size(b{nn,1},1),1,i)
% plot(b{1,1}(i,:))
% xlim ([0 10000])
% end
% end
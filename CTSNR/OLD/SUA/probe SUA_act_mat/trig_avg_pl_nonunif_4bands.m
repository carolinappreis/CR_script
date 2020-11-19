clear all
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/probe SUA_act_mat')
load ('NEW_SNR_cycle.mat','units_match','ecogbf_match') %%% ecog is filtered i differetn band of beta 1 ou of 4 to which spikes have locked better and were not uniformily distributed

units_match1=units_match(~cellfun('isempty',units_match));
ecogbf_match1=ecogbf_match(any(ecogbf_match,2),:);

ctx=ecogbf_match1;
[env_var]=burst_var(ctx);

data=cell2mat(units_match1);
Ecogfiltered=[];
for i=1:size(units_match1,1)  
       Ecogfiltered = [Ecogfiltered ; repmat(ecogbf_match1(i,:),size(units_match1{i,1},1),1)];
end

clearvars -except data Ecogfiltered 

srn=1000;


tt=0;
for j=1:size(data,1)
    
    clearvars -except j Ecogfiltered data srn rec_pa rec_npa count_b tt
    
    
        env=abs(hilbert(Ecogfiltered(j,:)));
        
        onset1=bursts(env);
        onset1=horzcat(onset1{:});
        onset=bursts_aligned(env,Ecogfiltered(j,:));
        onset=horzcat(onset{:});
        
        data_g=smoothdata(data(j,:),'gaussian',25);
        
        for jj=1:size(onset,2)
            if onset(jj)>200 && onset(jj)+200<length(data_g) && onset1(jj)>200 && onset1(jj)+200<length(data_g) 
                output_pa(jj,:)=data_g(onset(jj)-200:onset(jj)+200);
                output_npa(jj,:)= data_g(onset1(jj)-200:onset1(jj)+200);
            end
        end
    
        if max(zscore(sum(output_pa,1)))>1.96
            tt=tt+1;
            count_b(1,tt)=length(onset);
            rec_pa(tt,:)=sum(output_pa,1);
            rec_npa(tt,:)=sum(output_npa,1);
        end
end


region_pl=zscore(mean(rec_pa,1));
region_spl=zscore(std(rec_pa)./sqrt(size(rec_pa,1)));
region_npl=zscore(mean(rec_npa,1));
region_snpl=zscore(std(rec_npa)./sqrt(size(rec_npa,1)));


%%% just plot
time2=[1:401];
 color_s=[0.5 0.5 0.5];
 color_b=[0 0 0.5]; %BZ
%   color_b=[0.5 0 0]; %bz

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





%%%%%%%% stats and images on filterd and evelope rtiggered average before
%%%%%%%% and after beta burst onset


% % clear j
% % for j=1:size(rec_pa,1)
% %     filt_pa2(j,:)=filtfilt(b,a,rec_pa(j,:));
% %     filt_npa2(j,:)=filtfilt(b,a,rec_npa(j,:));
% %     env_pa2(j,:)=abs(hilbert(filt_pa2(j,:)));
% %     env_npa2(j,:)=abs(hilbert(filt_npa2(j,:)));
% %     stat_pa2(j,1:2)=[sum(env_pa2(j,1:200)) sum(env_pa2(j,201:400))];
% %     stat_npa2(j,1:2)=[sum(env_npa2(j,1:200)) sum(env_npa2(j,201:400))];
% % end

% % 
% % [p(1) h(1)]=ttest(stat_pa2(:,1),stat_pa2(:,2));
% % [p(2) h(2)]=ttest(stat_npa2(:,1),stat_npa2(:,2));
% % [p(3) h(3)]=ttest(stat_pa2(:,1),stat_npa2(:,1));
% % [p(4) h(4)]=ttest(stat_pa2(:,2),stat_npa2(:,2));

% % 
% % figure;
% % subplot(1,2,1)
% % plot(env_pa2')
% % hold on
% % plot(mean(env_pa2),'k','LineWidth',2)
% % xlim([0 400])
% % box('off')
% % title('phase-aligned envelope')
% % 
% % 
% % subplot(1,2,2)
% % plot(env_npa2')
% % hold on
% % plot(mean(env_npa2),'k','LineWidth',2)
% % xlim([0 400])
% % box('off')
% % title('non-phase-aligned envelope')
% % 
% % 
% % fig=figure;
% % subplot(1,2,1)
% % bar(1:2,[median(stat_pa2(:,1)) median(stat_pa2(:,2))],'FaceColor',[0.5 0.5 0.5],'FaceAlpha',.7,'EdgeColor','none','BarWidth',0.5);
% % hold on
% % plot(1,stat_pa2(:,1),'.')
% % plot(2,stat_pa2(:,2),'.')
% % box('off')
% % xlim([0.5 2.5])
% % ylim([0 60])
% % title('phase-aligned envelope')
% % txt=(sprintf('p= %d',( round(h(1),4))));
% % text(1,58,txt)
% % set(gca,'FontSize',12)
% % 
% % subplot(1,2,2)
% % 
% % bar(1:2,[median(stat_npa2(:,1)) median(stat_npa2(:,2))],'FaceColor',[0.5 0.5 0.5],'FaceAlpha',.7,'EdgeColor','none','BarWidth',0.5);
% % hold on
% % plot(1,stat_npa2(:,1),'.')
% % plot(2,stat_npa2(:,2),'.')
% % box('off')
% % xlim([0.5 2.5])
% % ylim([0 60])
% % title('non-phase-aligned envelope')
% % txt=(sprintf('p= %d',( round(h(2),4))));
% % text(1,58,txt)
% % set(gca,'FontSize',12)
% % 
% % fig.Units = 'centimeters';
% % fig.OuterPosition= [10, 10, 18, 15];
% % fig.Color='w';



% % [b,a]=butter(2,[20/(0.5*srn) 30/(0.5*srn)],'bandpass');
% % 
% % fig=figure;
% % subplot(1,2,1)
% % y2=filtfilt(b,a,region_pl); y1=y2+filtfilt(b,a,region_spl); y3=y2-filtfilt(b,a,region_spl);
% % y5=filtfilt(b,a,region_npl); y4=y5+filtfilt(b,a,region_snpl); y6=y5-filtfilt(b,a,region_spl);
% % p1=plot(time2, y2, 'LineWidth',1.5,'Color',color_b)
% % patch([time2 fliplr(time2)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.1],'EdgeColor','none')
% % patch([time2 fliplr(time2)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.1],'EdgeColor','none')
% % xline(200,'--',{'burst onset'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','Color',[0.5 0.5 0.5],'LineWidth',2)
% % xlim ([0 400])
% % ylim ([-5 5])
% % xticks([0:100:400])
% % xticklabels ({'-200','-100','0','100','200'})
% % legend([p1],{'phase-aligned'},'FontSize',12,'box','off','Location','northeast' )
% % box ('off')
% % xlabel ('Time (msec)')
% % ylabel('Firing-rate(z-score)')
% % set(gca,'FontSize',12)
% % 
% % 
% % subplot(1,2,2)
% % p2=plot(time2, y5, 'LineWidth',1.5,'Color',color_s)
% % patch([time2 fliplr(time2)], [y4 fliplr(y5)],[color_s],'FaceAlpha',[0.1],'EdgeColor','none')
% % patch([time2 fliplr(time2)], [y5 fliplr(y6)],[color_s],'FaceAlpha',[0.1],'EdgeColor','none')
% % xline(200,'--',{'burst onset'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','Color',[0.5 0.5 0.5],'LineWidth',2)
% % xlim ([0 400])
% % ylim ([-5 5])
% % xticks([0:100:400])
% % xticklabels ({'-200','-100','0','100','200'})
% % legend([p2],{'non-phase aligned'},'FontSize',12,'box','off','Location','northeast' )
% % box ('off')
% % xlabel ('Time (msec)')
% % ylabel('Firing-rate(z-score)')
% % set(gca,'FontSize',12)
% % 
% % fig.Units = 'centimeters';
% % fig.OuterPosition= [10, 10, 15, 10];
% % fig.Color='w';
% % 
% % 

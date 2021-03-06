clear all
%  cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\probe SUA_act_mat')
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/probe SUA_act_mat')
load('data_SUA_BZ.mat')

srn=1000;

% % % check firing rate units
% data_a=cell2mat(data_region);
% % % for ii=1:size(data_a,1)
% % %     data=data_a(ii,:);
% % %     
% % %     spkrate_1=[];
% % %     for i =1:srn:(length(data)-srn);
% % %         spkrate_1=[spkrate_1 numel(find(data(i:i+srn)==1))];
% % %     end
% % %     srate(ii,1)=mean(spkrate_1);
% % %     clear data
% % % end


[b,a]=butter(2,[5/(0.5*srn) 15/(0.5*srn)],'bandpass');
[bb,aa]=butter(2,[15/(0.5*srn) 35/(0.5*srn)],'bandpass');
for j=1:size(data_region,1)
    
    data=data_region{j,1};
    ecog=Ecog_region(j,:);
    clearvars -except bb aa j ii BZ b a stat_pa2 stat_npa2 env_pa2 env_npa2 filt_pa2 filt_npa2 Ecog_region data ecog data_region srn rec_pa1 rec_npa1 rec_pa2 rec_npa2 region_pl region_spl region_npl region_nspl subj color_b
    for ii=1:size(data,1)
        Ecogfiltered=filtfilt(bb,aa,ecog);
        env=abs(hilbert(Ecogfiltered));
        onset1=bursts(env);
        onset1=horzcat(onset1{:});
        onset=bursts_aligned(env,Ecogfiltered);
        onset=horzcat(onset{:});
        
        data_g=smoothdata(data(ii,:),'gaussian',25);
        
        for jj=1:size(onset,2)
            if onset(jj)>200 && onset(jj)+200<length(data_g)
                output_pa(ii,jj,:)=data_g(onset(jj)-200:onset(jj)+200);
                output_npa(ii,jj,:)= data_g(onset1(jj)-200:onset1(jj)+200);
                output_count(ii,jj,:)= data(onset1(jj)-200:onset1(jj)+200);
            end
        end
    end
        
    rec_pa1{j,1}=reshape(sum(output_pa,2),(size(sum(output_pa,2),1)),(size(sum(output_pa,2),3)));
    rec_npa1{j,1}=reshape(sum(output_npa,2),(size(sum(output_npa,2),1)),(size(sum(output_npa,2),3)));
    count{j,1}=reshape(sum(output_count,2),(size(sum(output_count,2),1)),(size(sum(output_count,2),3)));
    clear output_pa output_npa output_count
    
    rec_pa2(j,:)=mean(rec_pa1{j,1},1);
    rec_npa2(j,:)=mean(rec_npa1{j,1},1);
    filt_pa2(j,:)=filtfilt(b,a,(mean(rec_pa1{j,1},1)));
    filt_npa2(j,:)=filtfilt(b,a,(mean(rec_npa1{j,1},1)));
    env_pa2(j,:)=abs(hilbert(filt_pa2(j,:)));
    env_npa2(j,:)=abs(hilbert(filt_npa2(j,:)));
    stat_pa2(j,1:2)=[sum(env_pa2(1:200)) sum(env_pa2(201:400))];
    stat_npa2(j,1:2)=[sum(env_npa2(1:200)) sum(env_npa2(201:400))];
    clear data ecog
end
ttest(stat_pa2(:,1),stat_pa2(:,2))
ttest(stat_npa2(:,1),stat_npa2(:,2))

ttest(stat_pa2(:,1),stat_npa2(:,1))
ttest(stat_pa2(:,2),stat_npa2(:,2))

for i =1:size(rec_pa2,1)
    subplot(size(rec_pa2,1),1,i)
    plot(rec_pa2(i,:))
end

%%%RAW
region_pl=zscore(mean(rec_pa2,1));
region_spl=zscore(std(rec_pa2)./sqrt(size(rec_pa2,1)));
region_npl=zscore(mean(rec_npa2,1));
region_snpl=zscore(std(rec_npa2)./sqrt(size(rec_npa2,1)));

% %%FILTERED
% region_pl=zscore(mean(filt_pa2,1));
% region_spl=zscore(std(filt_pa2)./sqrt(size(filt_pa2,1)));
% region_npl=zscore(mean(filt_npa2,1));
% region_snpl=zscore(std(filt_npa2)./sqrt(size(filt_npa2,1)));

%%%ENVELOPE
% region_pl=zscore(mean(env_pa2,1));
% region_spl=zscore(std(env_pa2)./sqrt(size(filt_pa2,1)));
% region_npl=zscore(mean(env_npa2,1));
% region_snpl=zscore(std(env_npa2)./sqrt(size(filt_npa2,1)));


time2=[1:401];
color_s=[0.5 0.5 0.5];
% color_b=[0 0 0.5]; %BZ
   color_b=[0.5 0 0]; %BZ

fig=figure;
subplot(1,2,1)
y2=region_pl; y1=y2+region_spl; y3=y2-region_spl;
y5=region_npl; y4=y5+region_snpl; y6=y5-region_snpl;
p1=plot(time2, y2, 'LineWidth',1.5,'Color',color_b)
patch([time2 fliplr(time2)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.1],'EdgeColor','none')
patch([time2 fliplr(time2)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.1],'EdgeColor','none')
hold on
xline(200,'--',{'burst onset'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','Color',[0.5 0.5 0.5],'LineWidth',2)
xlim ([0 400])
ylim ([-5 5])
xticks([0:100:400])
xticklabels ({'-200','-100','0','100','200'})
% legend([p1],{'phase-aligned'},'FontSize',12,'box','off')
box ('off')
xlabel ('Time (msec)')
ylabel('Firing-rate(z-score)')
set(gca,'FontSize',12)

% title('Spike triggered average')

subplot(1,2,2)
p2=plot(time2, y5, 'LineWidth',1.5,'Color',color_s)
patch([time2 fliplr(time2)], [y4 fliplr(y5)],[color_s],'FaceAlpha',[0.1],'EdgeColor','none')
patch([time2 fliplr(time2)], [y5 fliplr(y6)],[color_s],'FaceAlpha',[0.1],'EdgeColor','none')
hold on
xline(200,'--',{'burst onset'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','Color',[0.5 0.5 0.5],'LineWidth',2)
xlim ([0 400])
ylim ([-5 5])
xticks([0:100:400])
xticklabels ({'-200','-100','0','100','200'})
legend([p2],{'non-phase aligned'},'FontSize',12,'box','off')
box ('off')
xlabel ('Time (msec)')
ylabel('Firing-rate(z-score)')

fig.Units = 'centimeters';
fig.InnerPosition= [10, 10, 16, 9];
fig.Color='w';
set(gca,'FontSize',12)


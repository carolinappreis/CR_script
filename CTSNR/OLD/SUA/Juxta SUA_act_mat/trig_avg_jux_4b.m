clear
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/Juxta SUA_act_mat/mat')
load('BZ_cycle')

data=units_match;
Ecogfiltered =ecogbf_match;

% ctx=Ecogfiltered;
% [env_var]=burst_var(ctx);

clearvars -except data Ecogfiltered 

tt=0;
for j=1:size(data,1)
    
    clearvars -except j BZ Ecogfiltered data srn rec_pa rec_npa count_b tt
    
    
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
%  color_b=[0 0 0.5]; %BZ
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
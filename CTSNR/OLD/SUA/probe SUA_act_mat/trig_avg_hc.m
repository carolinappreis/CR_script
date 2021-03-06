clear all
 cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\probe SUA_act_mat')
% cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/probe SUA_act_mat')
load('NEW_SNR_cycle.mat')

data_region=units_match(~cellfun('isempty',units_match));
Ecog_region=ecogbf_match(any(ecogbf_match,2),:);

srn=1000;

for j=1:size(data_region,1)

    clearvars -except j ii SNR Ecog_region  ecog data_region srn rec_pa1 rec_npa1 rec_pa2 rec_npa2 region_pl region_spl region_npl region_nspl subj color_b
    
    
    data=data_region{j,1};
    Ecogfiltered=Ecog_region(j,:);
    
    for ii=1:size(data,1)
        env=abs(hilbert(Ecogfiltered));
        onset1=bursts(env);
        onset1=horzcat(onset1{:});
        onset=bursts_aligned(env,Ecogfiltered);
        onset=horzcat(onset{:});
        
        data_g=smoothdata(data(ii,:),'gaussian',25);
        
        epoch=500;
        for jj=1:size(onset,2)
            if onset(jj)>epoch && onset(jj)+epoch<length(data_g)
                output_pa(ii,jj,:)=data_g(onset(jj)-epoch:onset(jj)+epoch);
            end
            if onset1(jj)>epoch && onset1(jj)+epoch<length(data_g)
                output_npa(ii,jj,:)= data_g(onset1(jj)-epoch:onset1(jj)+epoch);
                output_count(ii,jj,:)= data(onset1(jj)-epoch:onset1(jj)+epoch);
            end
        end
        rec_pa1{j,1}=reshape(sum(output_pa,2),(size(sum(output_pa,2),1)),(size(sum(output_pa,2),3)));
        rec_npa1{j,1}=reshape(sum(output_npa,2),(size(sum(output_npa,2),1)),(size(sum(output_npa,2),3)));
    end
    rec_pa2(j,:)=mean(rec_pa1{j,1},1);
    rec_npa2(j,:)=mean(rec_npa1{j,1},1);
    clear data ecog
end


region_pl=zscore(mean(rec_pa2,1));
region_spl=zscore(std(rec_pa2)./sqrt(size(rec_pa2,1)));
region_npl=zscore(mean(rec_npa2,1));
region_snpl=zscore(std(rec_npa2)./sqrt(size(rec_npa2,1)));



time2=[1:1001];
color_s=[0.5 0.5 0.5];
color_b=[0 0 0.5]; %SNR
% color_b=[0.5 0 0]; %bz

fig=figure;
subplot(1,2,1)
y2=region_pl; y1=y2+region_spl; y3=y2-region_spl;
y5=region_npl; y4=y5+region_snpl; y6=y5-region_snpl;
p1=plot(time2, y2, 'LineWidth',1.5,'Color',color_b)
patch([time2 fliplr(time2)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.1],'EdgeColor','none')
patch([time2 fliplr(time2)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.1],'EdgeColor','none')
% xlim ([0 400])
% ylim ([-5 5])
% xticks([0:100:400])
% xticklabels ({'-epoch','-100','0','100','epoch'})
xlim ([300 700])
ylim ([-5 5])
% xticks([0:100:400])
xticklabels ({'-200','-100','0','100','200'})
legend([p1],{'phase-aligned'},'FontSize',12,'box','off')
box ('off')
xlabel ('Time (msec)')
ylabel('Firing-rate(z-score)')


subplot(1,2,2)
p2=plot(time2, y5, 'LineWidth',1.5,'Color',color_s)
patch([time2 fliplr(time2)], [y4 fliplr(y5)],[color_s],'FaceAlpha',[0.1],'EdgeColor','none')
patch([time2 fliplr(time2)], [y5 fliplr(y6)],[color_s],'FaceAlpha',[0.1],'EdgeColor','none')
xlim ([300 700])
ylim ([-5 5])
% xticks([0:100:400])
xticklabels ({'-200','-100','0','100','200'})
legend([p2],{'non-phase aligned'},'FontSize',12,'box','off')
box ('off')
xlabel ('Time (msec)')
ylabel('Firing-rate(z-score)')

fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 15, 10];
fig.Color='w';

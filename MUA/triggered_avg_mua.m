clear all
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\MUA')
load('BZ_spikerate.mat') ; 

Ecog=BZ.ctx(BZ.beta_rats,:);
for i =1:size(BZ.beta_rats,1)
    for ii=1:size(BZ.beta_idx{BZ.beta_rats(i),1},2)
        hp=BZ.beta_idx{BZ.beta_rats(i),1}(1,ii);
        data_all{i,1}(ii,:)=BZ.sua{BZ.beta_rats(i),1}(hp,:);
        clear hp
    end
end


srn=1000;

for j=1:size(data_all,1)
    
    data=data_all{j,1};
    ecog=Ecog(j,:);
    clearvars -except j ii srn Ecog data ecog data_all  rec_pa1 rec_npa1 region_pl region_spl region_npl region_nspl color_b
    
    for ii=1:size(data,1)
        
        [b,a]=butter(2,[15/(0.5*srn) 35/(0.5*srn)],'bandpass');
        Ecogfiltered=filtfilt(b,a,ecog);
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
        rec_pa1{j,1}=reshape(sum(output_pa,2),(size(sum(output_pa,2),1)),(size(sum(output_pa,2),3)));
        rec_npa1{j,1}=reshape(sum(output_npa,2),(size(sum(output_npa,2),1)),(size(sum(output_npa,2),3)));
    end
    
    clear data
end
rec_pa=cell2mat(rec_pa1);
rec_npa=cell2mat(rec_npa1);

region_pl=zscore(mean(rec_pa,1));
region_spl=zscore(std(rec_pa)./sqrt(size(rec_pa,1)));
region_npl=zscore(mean(rec_npa,1));
region_snpl=zscore(std(rec_npa)./sqrt(size(rec_npa,1)));

time2=[1:401];
color_b=[0.2 0.5 0.5];
color_b=[0.5 0 0.5];
color_s=[0.5 0.5 0.5];



subplot(1,2,1)
y2=region_pl; y1=y2+region_spl; y3=y2-region_spl;
y5=region_npl; y4=y5+region_snpl; y6=y5-region_snpl;
p1=plot(time2, y2, 'LineWidth',1.5,'Color',color_b)
patch([time2 fliplr(time2)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.1],'EdgeColor','none')
patch([time2 fliplr(time2)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.1],'EdgeColor','none')
xlim ([0 400])
ylim ([-5 5])
xticks([0:100:400])
xticklabels ({'-200','-100','0','100','200'})
legend([p1],{'phase-aligned'},'FontSize',12)
box ('off')
%         xlabel ('Time (msec)')
%         ylabel('Firing-rate(z-score)')


subplot(1,2,2)
p2=plot(time2, y5, 'LineWidth',1.5,'Color',color_s)
patch([time2 fliplr(time2)], [y4 fliplr(y5)],[color_s],'FaceAlpha',[0.1],'EdgeColor','none')
patch([time2 fliplr(time2)], [y5 fliplr(y6)],[color_s],'FaceAlpha',[0.1],'EdgeColor','none')
xlim ([0 400])
ylim ([-5 5])
xticks([0:100:400])
xticklabels ({'-200','-100','0','100','200'})
legend([p2],{'non-phase aligned'},'FontSize',12)
box ('off')


% plot(rec_pa1{10,:}(2,:)') % very good cz sua
% plot(rec_pa1{1,:}')% snr example
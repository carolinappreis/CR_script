clear all
cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A2_Thal/MAT')
load('BZ_spikerate.mat') ;


% subj= BZ.beta_rats(ismember(BZ.beta_rats,SNR.beta_rats));
ii=1;
srn=1000;output_pa=cell(1,1);output_npa=cell(1,1);
for j=1:size(BZ.beta_rats,1)
    
    data= BZ.sua{BZ.beta_rats(j),1};
    
    for ii=1:size(data,1)
        
            clearvars -except output_pa BZ data srn output_count rec_pa rec_npa A all_regions ii j region_pl region_npl region_spl region_snpl output_pa output_npa ISI_all spikerate_all

        
        % cd('C:\Users\creis\Documents\GitHub\CR_script\A4_Thal\code')
        %         cd('/Users/Carolina/Documents/GitHub/CR_script/codes_thal/A4_Thal/code')
        [b,a]=butter(2,[15/(0.5*srn) 35/(0.5*srn)],'bandpass');
        Ecogfiltered=filtfilt(b,a,BZ.ctx(BZ.beta_rats(j),:));
        env=abs(hilbert(Ecogfiltered));
        onset1=bursts(env);
        onset1=horzcat(onset1{:});
        onset=bursts_aligned(env,Ecogfiltered);
        onset=horzcat(onset{:});
        
        data_g=smoothdata(data(ii),'gaussian',25);
        
        for nb=1:size(onset1,1)
            for jj=1:size(onset,2)
                if onset(jj)>200 && onset(jj)+200<length(data_g)
                    output_pa{j,:}(jj,:)=data_g(onset(jj)-200:onset(jj)+200);
                    output_npa{j,:}(jj,:)= data_g(onset1(jj)-200:onset1(jj)+200);
                    output_count{j,:}(jj,:)= data(onset1(jj)-200:onset1(jj)+200);
                end
            end
            rec_pa(j,:)=sum(output_pa{j,:},1);
            rec_npa(j,:)=sum(output_npa{j,:},1);
            
            region_pl{ii,1}=zscore(mean(rec_pa,1));
            region_spl{ii,1}=zscore(std(rec_pa)./sqrt(size(rec_pa,1)));
            region_npl{ii,1}=zscore(mean(rec_npa,1));
            region_snpl{ii,1}=zscore(std(rec_npa)./sqrt(size(rec_npa,1)));
        end
    end
end

time2=[1:401];
color_b=[0.2 0.5 0.5];
color_b=[0.5 0 0.5];
color_s=[0.5 0.5 0.5];
titles={'CZ','BZ'};
for b=1
    for j =2
        %         size(region_pl,1)
        subplot(1,2,1)
        y2=region_pl{j,b}; y1=y2+region_spl{j,b}; y3=y2-region_spl{j,b};
        y5=region_npl{j,b}; y4=y5+region_snpl{j,b}; y6=y5-region_snpl{j,b};
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
        
        title (titles(j),'FontSize',12)
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
    end
end


clear all
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
load ('data_all' , 'freq')

cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/Juxta SUA_act_mat/mat')
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\Juxta SUA_act_mat\mat')
load 'animal_region.mat' % comes from code pre_filegen_SUA_act

all_regions={[thal_VL];[thal_VA thal_VM]};
for ii=1:size(all_regions,1)
    for  j=1:size(all_regions{ii,:},2)
        
        clearvars -except output_count rec_pa rec_npa A all_regions ii j region_pl region_npl region_spl region_snpl output_pa output_npa ISI_all spikerate_all
        name=A(all_regions{ii,:}(j),1:(find(A((all_regions{ii,:}(j)),:)=='.')-1));
%         cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\Juxta SUA_act_mat\mat')
               cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/Juxta SUA_act_mat/mat')
        
        load(name);
        B=who;
        ctxchan=[];
        for rr=1:size(B,1)
            if ~isempty(min(find(B{rr}=='E')) & min(find(B{rr}=='E')) & min(find(B{rr}=='G')))
                ctxchan=rr;
            end
        end
        
        %--------------------------------------------------------------
        eval(['samprateold=1/' B{ctxchan} '.interval;']);
        eval(['WaveData(1,:)=' B{ctxchan} '.values;']);
        WaveData=double(WaveData);
        ts=timeseries(WaveData(1,:),0:(1/samprateold):((size(WaveData,2)-1)/samprateold));
        ts1=resample(ts,0:0.001:((size(WaveData,2)-1)/samprateold),'linear');
        WaveData_DC(1,:)=ts1.data;
        
        %--------------------------------------------------------------
        sr=1/unite.interval;
        srn=1000;
        dataold=unite.values';
        dataold=full(dataold);
        data=zeros(1,100000);
        timeold=0:1/sr:(size(dataold,2)-1)/sr;
        time=0:1/srn:(size(data,2)-1)/srn;
        spk_t=timeold(find(dataold==1));
        spk_tround=round(spk_t,3);
        nn=[];
        for i=1:length(spk_t)
            [ d, ix ] = min( abs( time-spk_tround(i) ) );
            nn=[nn ix];
        end
        data(nn)=1;
        %         data_all{ii,j}=data;
        data_ones=find(data==1);
        %---------------------------
        
        [b,a]=butter(2,[15/(0.5*srn) 35/(0.5*srn)],'bandpass');
        Ecogfiltered=filtfilt(b,a,WaveData_DC);
        env=abs(hilbert(Ecogfiltered));
        spkrate_1=[];
        for i =1:srn:(length(data)-srn);
            spkrate_1=[spkrate_1 numel(find(data(i:i+srn)==1))];
            ISI_all{ii,1}(1,j)=mean(diff(time(data==1)));
            %                 psi_all{ii,1}(j,:)=abs(mean(find(ang(data_ones)));
        end
        
        spikerate_all{ii,1}(1,j)=mean(nonzeros(spkrate_1));
        
        %-----------------------------
        onset1=bursts(env);
        onset1=horzcat(onset1{:});
        onset=bursts_aligned(env,Ecogfiltered);
        onset=horzcat(onset{:});
        %         onset=bursts_aligned(env,Ecogfiltered);
        %
        %         plot(time,Ecogfiltered)
        %         hold on
        %         plot(time(onset1{1,1}),env(onset1{1,1}),'bo','MarkerSize', 5)
        %         plot(time(onset{1,1}),env(onset{1,1}),'r.','MarkerSize', 10)
        %         plot(time,tt)
        %         plot(time,env)
        %-------------------------------
        
        data_g=smoothdata(data,'gaussian',25);
        
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
        %         for nb=1:size(onset1,1)
        %             for jj=1:size(onset{nb,1},2)
        %                 if onset{nb,1}(jj)>200 && onset{nb,1}(jj)+200<length(data_g)
        %                     output_pa{j,:}(jj,:)=data_g(onset{nb,1}(jj)-200:onset{nb,1}(jj)+200);
        %                     output_npa{j,:}(jj,:)= data_g(onset1{nb,1}(jj)-200:onset1{nb,1}(jj)+200);
        %                 end
        %             end
        %             rec_pa{nb,1}(j,:)=sum(output_pa{j,:},1);
        %             rec_npa{nb,1}(j,:)=sum(output_npa{j,:},1);
        %
        %             region_pl{nb,ii}=zscore(mean(rec_pa{nb,1},1));
        %             region_spl{nb,ii}=zscore(std(rec_pa{nb,1})./sqrt(size(rec_pa{nb,1},1)));
        %             region_npl{nb,ii}=zscore(mean(rec_npa{nb,1},1));
        %             region_snpl{nb,ii}=zscore(std(rec_npa{nb,1})./sqrt(size(rec_npa{nb,1},1)));
        %         end
    end
end

time2=[1:401];
color_b=[0.2 0.5 0.5];
color_b=[0.5 0 0];
color_s=[0.5 0.5 0.5];
titles={'CZ','BZ'};
for b=1
    for i =2
        fig=figure(1)
        subplot(1,2,1)
        y2=region_pl{i,b}; y1=y2+region_spl{i,b}; y3=y2-region_spl{i,b};
        y5=region_npl{i,b}; y4=y5+region_snpl{i,b}; y6=y5-region_snpl{i,b};
        p1=plot(time2, y2, 'LineWidth',1.5,'Color',color_b)
        patch([time2 fliplr(time2)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.1],'EdgeColor','none')
        patch([time2 fliplr(time2)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.1],'EdgeColor','none')
        xlim ([0 400])
        ylim ([-5 5])
        xticks([0:100:400])
        xticklabels ({'-200','-100','0','100','200'})
        legend([p1],{'phase-aligned'},'FontSize',12,'box','off')
        box ('off')
        xlabel ('Time (msec)')
        ylabel('Firing-rate (z-score)')
        
%         title (titles(i),'FontSize',12)
        subplot(1,2,2)
        p2=plot(time2, y5, 'LineWidth',1.5,'Color',color_s)
        patch([time2 fliplr(time2)], [y4 fliplr(y5)],[color_s],'FaceAlpha',[0.1],'EdgeColor','none')
        patch([time2 fliplr(time2)], [y5 fliplr(y6)],[color_s],'FaceAlpha',[0.1],'EdgeColor','none')
        xlim ([0 400])
        ylim ([-5 5])
        xticks([0:100:400])
        xticklabels ({'-200','-100','0','100','200'})
        legend([p2],{'non-phase aligned'},'FontSize',12,'box','off')
        box ('off') 
        xlabel ('Time (msec)')
ylabel('Firing-rate (z-score)')
    end
end

fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 15, 10];
fig.Color='w';

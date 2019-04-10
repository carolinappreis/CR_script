clear all
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
% cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
load ('data_all' , 'freq')

% cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/Juxta SUA_act_mat/mat')
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\Juxta SUA_act_mat\mat')
load 'animal_region.mat' % comes from code pre_filegen_SUA_act
% load 'non_repeat_animal_region.mat'
% load 'single_subj_VM_VL.mat'
% all_regions={thal_VL; thal_VM};
all_regions={thal_VA; thal_VL; thal_VM};
for ii=1:size(all_regions,1)
    for  j=1:size(all_regions{ii,:},2)
        
        clearvars -except A all_regions ii j region_pl region_npl region_spl region_snpl output_pa output_npa ISI_all spikerate_all
        name=A(all_regions{ii,:}(j),1:(find(A((all_regions{ii,:}(j)),:)=='.')-1));
        load(name);
        B=who;
        ctxchan=[];
        for rr=1:size(B,1)
            if ~isempty(min(find(B{rr}=='E')) & min(find(B{rr}=='E')) & min(find(B{rr}=='G')))
                ctxchan=rr;
            end
        end
        
        
        eval(['samprateold=1/' B{ctxchan} '.interval;']);
        eval(['WaveData(1,:)=' B{ctxchan} '.values;']);
        WaveData=double(WaveData);
        ts=timeseries(WaveData(1,:),0:(1/samprateold):((size(WaveData,2)-1)/samprateold));
        ts1=resample(ts,0:0.001:((size(WaveData,2)-1)/samprateold),'linear');
        WaveData_DC(1,:)=ts1.data;
        
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
        data_all{ii,j}=data;
        data_ones=find(data==1);
        %---------------------------
        
        
        [b,a]=butter(2,[15/(0.5*srn) 35/(0.5*srn)],'bandpass');
        Ecogfiltered=filtfilt(b,a,WaveData_DC);
        env=abs(hilbert(Ecogfiltered));
        threshold=prctile(env,75);
        tt(size(env,1):size(env,2))=threshold;
        indexexceed=find(env>threshold);
        diffindex=diff(indexexceed);
        pnts=find(diffindex>1);
        begin=indexexceed(pnts+1);
        ending=indexexceed(pnts);
        begin2=[indexexceed(1) begin];
        ending2=[ending indexexceed(end)];
        duration=ending2-begin2;
        mean_b=mean(duration);
        ind_b=[];
        for i=1:(length(begin2))
            if (ending2(i)-begin2(i))>=100 % min duration of bursts
                ind_b=[ind_b i];
            end
        end
        
        if ~isempty (ind_b)
            begin3=begin2(ind_b);
            ending3=ending2(ind_b);
            
            space_betb=200; % min space between bursts
            ind_b1=[];
            for i=1:(length(begin3)-2)
                if (begin3(i+1)-ending3(i))>=space_betb && (begin3(i+2)-ending3(i+1))>=space_betb
                    ind_b1=[ind_b1 i+1];
                end
            end
            if (begin3(2)-ending(1))>=space_betb
                ind_b1= [1 ind_b1];
            end
            if (begin3(length(begin3))-ending3(length(begin3)-1))>=space_betb
                ind_b1=[ind_b1 length(begin3)];
            end
            onset1=begin3(ind_b1);
            offset1=ending3(ind_b1);
            
            %-----
            spkrate_1=[];
            for i =1:srn:(length(data)-srn);
                spkrate_1=[spkrate_1 numel(find(data(i:i+srn)==1))];
                ISI_all{ii,1}(1,j)=mean(diff(time(data==1)));
                %                 psi_all{ii,1}(j,:)=abs(mean(find(ang(data_ones)));
            end
            
            %-----
            spikerate_all{ii,1}(1,j)=mean(nonzeros(spkrate_1));
            
            
            [maxvalM,maxidxM] = findpeaks(Ecogfiltered);
            
            pre_onset=cell(1,1);
            for b = 1:length(onset1)
                for p=1:length(maxidxM)
                    if min(abs(onset1(b)-maxidxM(p)))<=30;
                        pre_onset{1,b}=p;
                    end
                end
            end
            pre_onset=cell2mat(pre_onset);
            onset=maxidxM(pre_onset);
            
            plot(time,Ecogfiltered)
            hold on
            plot(time(onset1),env(onset1),'bo','MarkerSize', 5)
            plot(time(onset),env(onset),'r.','MarkerSize', 10)
            plot(time,tt)
            plot(time,env)
            %-------------------------------
            data_g=data;
            
            for i=1:size(data_ones,2);
                if data_ones(i)>15
                    x=time(data_ones(i))-0.015:0.001:time(data_ones(i))+0.015;
                    data_g(1,data_ones(i)-15:data_ones(i)+15)=gaussmf(x,[0.005 time(data_ones(i))]);
                end
%                     plot(x,data_g(1,data_ones(i)-15:data_ones(i)+15))
            end
            
            for jj=1:length(onset)
                if onset(jj)>200 && onset(jj)+200<length(data_g)
                    output_pa{j,:}(jj,:)= data_g(onset(jj)-200:onset(jj)+200);
                    %                 output_pa{j,:}(jj,:)= data(onset1(jj)-200:onset1(jj)+200);
                end
            end
            for jj=1:length(onset1)
                if onset1(jj)>200 && onset1(jj)+200<length(data_g)
                    output_npa{j,:}(jj,:)= data_g(onset1(jj)-200:onset1(jj)+200);
                end
            end
            
            for i=1:size(output_pa,1);
                rec_pa(i,:)=sum(output_pa{i,1},1);
                rec_npa(i,:)=sum(output_npa{i,1},1);
            end
        end
    end
        region_pl(ii,:)=mean(rec_pa,1);
        region_spl(ii,:)=std(rec_pa)./sqrt(size(rec_pa,1));
        region_npl(ii,:)=mean(rec_npa,1);
        region_snpl(ii,:)=std(rec_npa)./sqrt(size(rec_npa,1));
end
time2=[1:401];
color_b=[0.2 0.5 0.5];
color_s=[0 0 0.5];
titles={'VA','VL','VM'};
for i =1:size(region_pl,1)
    subplot(3,1,i)
    y2=region_pl(i,:); y1=y2+region_spl(i,:); y3=y2-region_spl(i,:);
    y5=region_spl(i,:); y4=y5+region_spl(i,:); y6=y5-region_spl(i,:);
    p1=plot(time2, y2, 'LineWidth',1.5,'Color',color_b)
    patch([time2 fliplr(time2)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
    patch([time2 fliplr(time2)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
    hold on
    p2=plot(time2, y5, 'LineWidth',1.5,'Color',color_s)
    patch([time2 fliplr(time2)], [y4 fliplr(y5)],[color_s],'FaceAlpha',[0.2],'EdgeColor','none')
    patch([time2 fliplr(time2)], [y5 fliplr(y6)],[color_s],'FaceAlpha',[0.2],'EdgeColor','none')
    xlim ([0 400])
    ylim ([0 25])
    xticks([0:100:400])
    xticklabels ({'-200','-100','0','100','200'})
    box ('off')
    title (titles(i),'FontSize',14)
end

legend([p1 p2],{'phase-aligned','non-phase aligned'},'FontSize',12)
box('off')
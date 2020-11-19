clear all
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\probe SUA_act_mat')
% cd('/Users/Carolina/Documents/GitHub/CRcode/codes_thal/SUA/probe SUA_act_mat')
 cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/probe SUA_act_mat')
load('data_SUA_BZ.mat')
srn=1000;
for ii=1:size(data_all,1)
    for  j=1:size(data_all{ii,:},1)
        clearvars -except A lesion ii j r rec_pl rec_npl rec_spl rec_snpl output_pa output_npa spikerate_all ISI_all data_all Ecog_all srn time
        data= data_all{ii,1}(j,:) ;
        data_ones=find(data_all{ii,1}(j,:)==1);
        [b,a]=butter(2,[15/(0.5*srn) 30/(0.5*srn)],'bandpass');
        Ecogfiltered=filtfilt(b,a,Ecog_all(ii,:));
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
            
            %---------
            spkrate_1=[];
            for i =1:srn:(length(data)-srn);
                spkrate_1=[spkrate_1 numel(find(data(i:i+srn)==1))];
                ISI_all{ii,1}(1,j)=mean(diff(time(data==1)));
                %                 psi_all{ii,1}(j,:)=abs(mean(find(ang(data_ones)));
            end
            spikerate_all{ii,1}(1,j)=mean(nonzeros(spkrate_1));
            %--------------
            %             mean(nonzeros(spkrate_1))-std(nonzeros(spkrate_1))>=13 && mean(nonzeros(spkrate_1))+std(nonzeros(spkrate_1))<35
            if mean(nonzeros(spkrate_1))-std(nonzeros(spkrate_1))>=13 && mean(nonzeros(spkrate_1))+std(nonzeros(spkrate_1))<35
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
%                 
%                             plot(time,Ecogfiltered)
%                             hold on
%                             plot(time(onset1),env(onset1),'bo','MarkerSize', 5)
%                             plot(time(onset),env(onset),'r.','MarkerSize', 10)
%                             plot(time,tt)
%                             plot(time,env)
                %-------------------------------
                data_g=data;
                
                for i=1:size(data_ones,2);
                    if data_ones(i)>15
                        x=time(data_ones(i))-0.015:0.001:time(data_ones(i))+0.015;
                        data_g(1,data_ones(i)-15:data_ones(i)+15)=gaussmf(x,[0.005 time(data_ones(i))]);
%                         plot(x,data_g(1,data_ones(i)-15:data_ones(i)+15))
                    end
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
                    cont_pa(i,:)=nansum(output_pa{i,1},1);
                    cont_npa(i,:)=nansum(output_npa{i,1},1);
                end
                rec_pl(ii,1:401)=nanmean(cont_pa,1);
                rec_spl(ii,1:401)=nanstd(cont_pa)./sqrt(size(cont_pa,1));
                rec_npl(ii,1:401)=nanmean(cont_npa,1);
                rec_snpl(ii,1:401)=nanstd(cont_npa)./sqrt(size(cont_npa,1));
            end
        end
    end
    %     if mean(spikerate_all{ii,1})-std(spikerate_all{ii,1})>=13 && mean(spikerate_all{ii,1})+std(spikerate_all{ii,1})<=35
    %         r(ii,1)=1;
    %     else
    %         r(ii,1)=NaN;
    %     end
    %
    
end


rec_pl1 = rec_pl(any(rec_pl,2),:);
rec_spl1 = rec_spl(any(rec_spl,2),:);
rec_npl1 = rec_npl(any(rec_npl,2),:);
rec_snpl1 = rec_snpl(any(rec_snpl,2),:);

in=[];
for i=1:size(rec_pl,1)
    if rec_pl(i,1)~=0
        
        in=[in i];
    end
end

% rats_id=str2num(A(lesion(in),4:6))
% unNum     = unique(rats_id(:))
% [n,bin]   = histc(rats_id(:),unNum)

% n=0;
% for i=20
% plot(time2,mean(rec_pl(5:10,:)))
% hold on
% plot(time2,mean(region_npl(5:10,:)))
% hold off
% end



time2=[1:401];
color_b=[0.2 0.5 0.5];
color_s=[0 0 0.5];
for i =1:size(rec_pl1,1)
    subplot(size(rec_pl1,1),1,i)
    y2=rec_pl1(i,:); y1=y2+rec_spl1(i,:); y3=y2-rec_spl1(i,:);
    y5=rec_npl1(i,:); y4=y5+rec_snpl1(i,:); y6=y5-rec_snpl1(i,:);
    p1=plot(time2, y2, 'LineWidth',1.5,'Color',color_b);
    patch([time2 fliplr(time2)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
    patch([time2 fliplr(time2)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
    hold on
    p2=plot(time2, y5, 'LineWidth',1.5,'Color',color_s);
    patch([time2 fliplr(time2)], [y4 fliplr(y5)],[color_s],'FaceAlpha',[0.2],'EdgeColor','none')
    patch([time2 fliplr(time2)], [y5 fliplr(y6)],[color_s],'FaceAlpha',[0.2],'EdgeColor','none')
    xlim ([0 400])
    % %     ylim ([0 10])
    xticks([0:100:400])
    xticklabels ({'-200','-100','0','100','200'})
    box ('off')
    % %         title (titles(i))
    hold off

end


    legend([p1 p2],{'phase-aligned','non-phase aligned'})
    box('off')


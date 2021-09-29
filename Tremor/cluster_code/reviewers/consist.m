clear; close
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cluster_out_mc.mat','out');
m_ax=1;% change if the main axis is not always 1 - replace by array of main axes
load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','stone','squash','sapphire','azure');
color_b1=[blushred; aegean; stone];
close all
clearvars -except m_ax color_b1 out
for iii=10
    %     1:size(out.start_c,1)
    %%% for type=1:2
% %     close all
    type=1; feature={'out.change_c';'out.fchange'};
    dr=eval([(num2str(feature{type,1})) '{' num2str(iii) ',2}{m_ax,1}']);
    
    data=nanmedian(dr);
    nostim=eval(squeeze([num2str(feature{type,1}) '{' num2str(iii) ',1}(m_ax,:)']));
    n=[];
    n=[n ; data(find(data > prctile(nostim, 99.7917) | data< prctile(nostim, 0.2083)))];
    n_phases=numel(find(data > prctile(nostim, 99.7917) | data< prctile(nostim, 0.2083)));
    if ~isempty(find(data > prctile(nostim, 99.7917) | data< prctile(nostim, 0.2083)))
        phase=find(data > prctile(nostim, 99.7917) | data< prctile(nostim, 0.2083));
    else
        phase=NaN;
    end
    
    if (~isnan(phase))
        for p=1:length(phase)
%              dat=out.env_seg{iii,phase(p)-1};
            dat=out.zenv_seg{iii,phase(p)};
% %             f1=figure(p)
            for gg=1:size(dat,1)
%                 subplot(2,5,gg)
                x=1:size(dat,2);
%                 y=smooth((dat(gg,:).*(10*9.81/0.5)));
                y=smooth(((dat(gg,:)-mean(dat(gg,1:1000)))./mean(dat(gg,1:1000))));
                data_all(gg,:)=y;
% % %                 plot(x,y,'k','LineWidth',1)
% % %                 hold on
% % %                 mdp1=fitlm(x,y,'poly2');
% % %                 plot(mdp1.Fitted,'r','LineWidth',1.5)
% % %                 n=mdp1.Rsquared.Adjusted;
% % %                 filename=['adjr^2=',num2str(round(n,2))];
% % %                 xlim([0 5000])
% % %                 xticks([0:1000:5000])
% % %                 xticklabels({'0','1','2','3','4','5'})
% % % %                                 ylim([0 16])
% % % %                                 yticks(0:4:16)
% % % %                                 ylim([0 5])
% % %                 set(gca,'FontSize',12)
% % %                 aa=get(gca,'ylim');
% % % %                 yticks([(aa(1)):(((aa(2)-aa(1))/2)):(aa(2))]);
% % %                 t = text(3500,aa(2),filename);
% % %                 s = t.FontSize;
% % %                 t.FontSize = 10;
% % %                 ylabel('tremor severity (m/s^2)')
% % %                 xlabel ('time (s)')
% % %                 box('off')  
            end
% % %             set(f1,'color','w');
% % %             f1.OuterPosition= [1,339,1440,539];
            
            
            
            f1=figure(5)
            subplot(1,4,p)
            y2=smooth(nanmean(data_all))';
            y1=(y2+smooth(((nanmean(data_all)+nanstd(data_all))./sqrt(size(data_all,1))))');
            y3=(y2-smooth(((nanmean(data_all)+nanstd(data_all))./sqrt(size(data_all,1))))');
            time=1:6000;
            patch([time fliplr(time)], [y1 fliplr(y2)],'r','FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
            hold on
            patch([time fliplr(time)], [y2 fliplr(y3)],'r','FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
            plot(time,y2,'r','LineWidth',2)
            ylim([-0.5 0.5])
            yticks(-0.5:0.25:0.5)
            xlim([0 6000])
            xticks([0:1000:6000])
            xticklabels({'-1','0','1','2','3','4','5'})
            set(gca,'FontSize',12)
            ylabel('Normalised tremor severity (a.u)')
            xlabel ('time (s)')
            box('off')
            xline(1000,'--',{'stim ON'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','LineWidth',2,'Color',[0.5 0.5 0.5])
            set(f1,'color','w');
            
        end
        f1.OuterPosition= [1,100,1000,300];
    end
    
    
end


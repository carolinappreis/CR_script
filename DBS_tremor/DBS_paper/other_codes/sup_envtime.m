function [f1]=sup_envtime(out,match_ax,spiral)
type=1;
for iii=1:size(out.start_c,1)
    
    feature={'out.mod_amp';'out.mod_freq'}; %%% non stim of freq has not been calculated, uncomment that bit to see frc
    dr=eval([(num2str(feature{type,1})) '{' num2str(iii) ',2}{match_ax(2,iii,1),1}']);
    data=nanmedian(dr);
    
    %%%% main axis only
    pp1=squeeze(out.seg_zenv{iii,1}{match_ax(2,iii,1),1}(:,(find(data==(min(data)))),:));
    dum=find(~isnan(pp1(:,1)));
    for i=1:length(dum)
         pp(i,:)=(pp1(i,:)-(mean(pp1(i,1:1000))))./(mean(pp1(i,1:1000)));
    end
    
    f1=figure(type+2)
    subplot(1,4,iii)
    y2=nanmean(pp);
    y1=y2+((nanmean(pp)+nanstd(pp))./sqrt(size(pp,1)));
    y3=y2-((nanmean(pp)+nanstd(pp))./sqrt(size(pp,1)));
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
    
    
    %%%% all trials per subject
    % f1=figure(type+2)
    % for i=1:length(find(~isnan((pp(:,1)))))
    %     subplot(2,5,i)
    %     plot(smooth(pp(i,:)),'k','LineWidth',2)
    %     xlim([0 6000])
    %     xticks([0:1000:6000])
    %     xticklabels({'-1','0','1','2','3','4','5'})
    %     ylim([0 4])
    %     yticks(0:1:4)
    %     set(gca,'FontSize',12)
    %     aa=get(gca,'ylim');
    %     ylabel('tremor severity (m/s^2)')
    %     xlabel ('time (s)')
    %     box('off')
    %     xline(1000,'--',{'stim ON'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','LineWidth',2,'Color',[0.5 0.5 0.5])
    % end
    %     f1.OuterPosition= [1,150,750,450];
    
    
    %%%  all axes
    
    cl=[[0 0 0];[1 0 0];[0 0 1]];
    for n=1:3
        feature={'out.mod_amp';'out.mod_freq'}; %%% non stim of freq has not been calculated, uncomment that bit to see frc
        dr=eval([(num2str(feature{type,1})) '{' num2str(iii) ',2}{match_ax(2,iii,n),1}']);
        data=nanmedian(dr);
        
        pp1=squeeze(out.seg_zenv{iii,1}{match_ax(2,iii,n),1}(:,(find(data==(min(data)))),:));
        dum=find(~isnan(pp1(:,1)));
        for i=1:length(dum)
            pp(i,:)=(pp1(i,:)-(mean(pp1(i,1:1000))))./(mean(pp1(i,1:1000)));
        end
        clear dum
        
        f2=figure(type+3)
        subplot(1,4,iii)
        y2=nanmean(pp);
        y1=y2+((nanmean(pp)+nanstd(pp))./sqrt(size(pp,1)));
        y3=y2-((nanmean(pp)+nanstd(pp))./sqrt(size(pp,1)));
        time=1:6000;
        patch([time fliplr(time)], [y1 fliplr(y2)],[cl(n,:)],'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off')
        hold on
        patch([time fliplr(time)], [y2 fliplr(y3)],[cl(n,:)],'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off')
        plot(time,y2,'Color',cl(n,:),'LineWidth',1)
        ylim([-0.5 0.5])
        yticks(-0.5:0.25:0.5)
        xlim([0 6000])
        xticks([0:1000:6000])
        xticklabels({'-1','0','1','2','3','4','5'})
        %     ylim([0 4])
        %     yticks(0:1:4)
        set(gca,'FontSize',12)
        ylabel('Normalised tremor severity (a.u)')
        xlabel ('time (s)')
        box('off')
        xline(1000,'--',{'stim ON'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','LineWidth',2,'Color',[0.5 0.5 0.5])
        set(f2,'color','w');
    end
%     if iii==4 && type==1 && spiral==0
%         [p1]=pt4(out,iii,type,match_ax);
%     end
end
f1.OuterPosition= [1,100,1000,300];
f2.OuterPosition= [1,100,1000,300];

end
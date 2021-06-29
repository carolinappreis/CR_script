function [p1]=plotting(out)


m_ax=1;% change if the main axis is not always 1 - replace by array of main axes
load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','stone','squash','sapphire','azure');
color_b1=[aegean;stone;blushred];


%%%  plot arc's main axis with significant threshold\
for iii=1:size(out.start_c,1)
    type=1;
    feature={'out.mod_amp';'out.mod_freq'}; %%% non stim of freq has not been calculated, uncomment that bit to see frc
    p=type;
    dr=eval([(num2str(feature{type,1})) '{' num2str(iii) ',2}{m_ax,1}']);
    data=nanmedian(dr);
    sup(iii,type)=min(data);
    
    nostim=eval(squeeze([num2str(feature{type,1}) '{' num2str(iii) ',1}(m_ax,:)']));
    n=[];
    n=[n ; data(find(data > prctile(nostim, 99.7917) | data< prctile(nostim, 0.2083)))];
    n_phases=numel(find(data > prctile(nostim, 99.7917) | data< prctile(nostim, 0.2083)));
    
    if ~isempty(find(data > prctile(nostim, 99.7917) | data< prctile(nostim, 0.2083)))
        phase=find(data > prctile(nostim, 99.7917) | data< prctile(nostim, 0.2083));
    else
        phase=NaN;
    end
    
    p1=figure(type)
    subplot(1,4,iii)
    y2=data;
    y1=prctile(dr,75);
    y3=prctile(dr,25);
    time=0:30:330;
    patch([time fliplr(time)], [y1 fliplr(y2)],[color_b1(p,:)],'FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
    hold on
    patch([time fliplr(time)], [y2 fliplr(y3)],[color_b1(p,:)],'FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
    stem(time,y2,'.', 'LineWidth',4,'MarkerSize',20,'Color',color_b1(p,:))
    hold on
    yline(0)
    
    plot(time,dr,'r.','MarkerSize',5)
    if ~isnan(phase)
        plot(time(phase),y2(phase),'.','Color',color_b1(3,:),'MarkerSize',25)
        %                 yline(prctile(nostim,99.7917),'k--','LineWidth',1)
        %                 yline(prctile(nostim,0.2083),'k--','LineWidth',1)
    end
    
    
    if type==2
        ylim([-0.5 0.5])
        yticks(-0.5:0.25:0.5)
        ylabel({'Change in tremor frequency (Hz)'})
    else
        %          ylim([-1.5 1.5])
        ylim([-1 1])
        ylabel({'Change in tremor severity'})
    end
    xlim([-5 335])
    xticks([0:30:330])
    box('off')
    xlabel({'Stimulation phase (degrees)'})
    set(gca,'FontSize',10)
    set(gca,'FontName','Arial','XTickLabelRotation',45)
    set(p1,'color','w');
    title(sprintf('patient %d',(iii)))  

    %%%  temporal evolution acc in most sup phase 

    pp=squeeze(out.seg_zenv{1,1}{1,1}(:,(find(data==(min(data)))),:));
    f1=figure(type+2)
    subplot(1,4,iii)
    plot(smooth(nanmedian(pp)),'r','LineWidth',2)
    xlim([0 6000])
    xticks([0:1000:6000])
    xticklabels({'-1','0','1','2','3','4','5'})
    ylim([0 4])
    yticks(0:1:4)
    set(gca,'FontSize',12)
    aa=get(gca,'ylim');
    ylabel('tremor severity (m/s^2)')
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
end

p1.OuterPosition= [1,100,1000,300];
f1.OuterPosition= [1,100,1000,300];
  
end
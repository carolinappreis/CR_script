function [nc]=dbs_plots_nc(out)
for iii=1:4
    
    load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','stone','squash');
    color_b1=[blushred; aegean; stone];
    
    
    m_ax=1;% change if the main axis is not always 1 - replace by array of main axes
    p=2;
    dr=squeeze(out.change_nc{iii,2}{m_ax,1});
    % nostim=squeeze(nostim1(iii,m_ax,:));
    nostim=squeeze(out.change_nc{iii,1}(m_ax,:));
    
    n=[];
    data=nanmedian(dr);
    n=[n ; data(find(data > prctile(nostim, 99.7917) | data< prctile(nostim, 0.2083)))];
    n_phases=numel(find(data > prctile(nostim, 99.7917) | data< prctile(nostim, 0.2083)));
    if ~isempty(find(data > prctile(nostim, 99.7917) | data< prctile(nostim, 0.2083)))
        phase=find(data > prctile(nostim, 99.7917) | data< prctile(nostim, 0.2083));
    else
        phase=NaN;
    end
    
    nc=figure(20)
    subplot(1,4,iii)
    y2=nanmedian(dr,1);
    y1=prctile(dr,75);
    y3=prctile(dr,25);
    time=0:30:330;
    patch([time fliplr(time)], [y1 fliplr(y2)],[color_b1(p,:)],'FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
    hold on
    patch([time fliplr(time)], [y2 fliplr(y3)],[color_b1(p,:)],'FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
    % plot(y2,'.', 'MarkerSize',20,'Color',color_b1(p,:))
    stem(time,y2,'.', 'LineWidth',1,'MarkerSize',10,'Color',color_b1(p,:))
    yline(0)
    if ~isnan(phase)
        plot(time(phase),y2(phase),'.','Color',color_b1(1,:),'MarkerSize',10)
    end
    yline(prctile(nostim,99.7917),'k--','LineWidth',1)
    yline(prctile(nostim,0.2083),'k--','LineWidth',1)
    ylim([-1 1])
    xlim([-5 335])
    xticks([0:30:330])
    box('off')
    
    ylabel({'Change in tremor severity'})
    xlabel({'Stimulation phase (degrees)'})
    set(gca,'FontSize',12)
    set(gca,'FontName','Arial','XTickLabelRotation',45)
    set(nc,'color','w');
end
end
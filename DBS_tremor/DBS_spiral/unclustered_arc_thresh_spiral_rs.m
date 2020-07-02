
load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/aux_out_spiral.mat')
m_ax=1;% change if the main axis is not always 1 - replace by array of main axes
load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','stone','squash','sapphire','azure');
color_b1=[blushred; aegean; stone];

% for iii=1:4
%     
%     rep = 10;
%     for ax=1:3
%         out.change_c{iii,2}{1,ax}=tt1_all{iii,ax};
%         bb=change_bl(iii,ax,:);
%         dum = bb(randi(length(bb), 1e6, rep)); clear bb
%         out.change_c{iii,1}(ax,:)=nanmedian(dum,2); clear dum
%         clear dum dum2 baseline3
%     end
%     
% end


for iii=1:4;
     p1=figure(1)
    type=1;
    feature={'out.change_c';'out.changef'};
    p=2;
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

   subplot(1,4,iii)
    
    y2=data;
    y1=prctile(dr,75);
    y3=prctile(dr,25);
    time=0:30:330;
    patch([time fliplr(time)], [y1 fliplr(y2)],[color_b1(p,:)],'FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
    hold on
    patch([time fliplr(time)], [y2 fliplr(y3)],[color_b1(p,:)],'FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
    % plot(y2,'.', 'MarkerSize',20,'Color',color_b1(p,:))
    stem(time,y2,'.', 'LineWidth',4,'MarkerSize',20,'Color',color_b1(p,:))
    yline(0)
    if ~isnan(phase)
        plot(time(phase),y2(phase),'.','Color',color_b1(1,:),'MarkerSize',25)
    end
%             yline(prctile(nostim,99.7917),'k--','LineWidth',1)
%             yline(prctile(nostim,0.2083),'k--','LineWidth',1)
    ylim([-1 1])
    xlim([-5 335])
    xticks([0:30:330])
    box('off')
    ylabel({'Change in tremor severity'})
    xlabel({'Stimulation phase (degrees)'})
    set(gca,'FontSize',9)
    set(gca,'FontName','Arial','XTickLabelRotation',45)
    set(p1,'color','w');
   
end
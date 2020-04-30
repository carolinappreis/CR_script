
clear all; close all;
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/data_matchcluster')
load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','squash');
cl=blushred;

n=[];
for i=1:size(tt1,1)
    %     f1 = figure(i)
    for axis = 1:size(tt1, 2)
        
        data=squeeze(tt1(i, axis,:));
        n=[n ; data(find(data > prctile(nostim(i, axis, :), 99.7917) | data< prctile(nostim(i, axis, :), 0.2083)))];
        n_phases(i,axis)=numel(find(data > prctile(nostim(i, axis, :), 99.7917) | data< prctile(nostim(i, axis, :), 0.2083)));
        if ~isempty(find(data > prctile(nostim(i, axis, :), 99.7917) | data< prctile(nostim(i, axis, :), 0.2083)))
            phase{i,axis}=find(data > prctile(nostim(i, axis, :), 99.7917) | data< prctile(nostim(i, axis, :), 0.2083));
        else
            phase{i,axis}=NaN;
        end
        
%         subplot(1, 3, axis)
%         bar(0:30:330, squeeze(tt1(i, axis,:)),'FaceColor',cl,'EdgeColor',cl)
%         hold on
%         yline(prctile(nostim(i, axis, :), 99.7917),'k--','LineWidth',1)
%         yline(prctile(nostim(i, axis, :), 0.2083),'k--','LineWidth',1)
%         ylim([-1 1])
%         ylabel('Change in tremor severity')
%         xlabel('Stimulation phase (degrees)')
%         set(gca,'XTickLabelRotation',45)
%         box('off')
        
    end
    
    f2= figure(20)
    load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','stone','squash');
    color_b1=[blushred; aegean; stone];
    p=2;
    %ARC sem/std patch
    subplot(2,5,i)
    dr=squeeze(tt_all{i,1}(:,:,main_clust(i)));
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
    if ~isnan(phase{i,1})
    plot(time(phase{i,1}),y2(phase{i,1}),'.','Color',color_b1(1,:),'MarkerSize',10)
    end
    ylim([-(max(y1)) max(y1)])
    xlim([-5 335])
    xticks([0:30:330])
    box('off')
    ylabel({'Change in tremor severity'})
    xlabel({'Stimulation phase (degrees)'})
    set(gca,'FontSize',12)
    set(gca,'FontName','Arial','XTickLabelRotation',45)
    
    
    %     f1.Units = 'centimeters';
    %     f1.OuterPosition= [10, 10, 30, 8];
    %     set(f1,'color','w');
    f2.Units = 'centimeters';
    f2.OuterPosition= [0, 0, 40, 15];
    set(f2,'color','w');
end



% sig_calc=size(tt1,1)*size(tt1,2)*size(tt1,3)*0.05;
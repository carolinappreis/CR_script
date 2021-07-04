function[avg_sup,f1]=three_by4(out,color_b1)
sup=[];cc=0;
type=2;
for iii=1:size(out.start_c,1)
    feature={'out.mod_amp';'out.mod_freq'}; %%% non stim of freq has not been calculated, uncomment that bit to see frc
    p=type;
    p1=figure(iii);
    for rr=1:3
        clearvars -except out match_ax color_b1 sup cc type iii rr feature p p1
        cc=cc+1;
        dr=eval([(num2str(feature{type,1})) '{' num2str(iii) ',2}{rr,1}']);
        data=nanmedian(dr);
        sup(type,iii,rr)=min(data);
        nostim=eval(squeeze([num2str(feature{type,1}) '{' num2str(iii) ',1}(rr,:)']));
        n=[];
        n=[n ; data(find(data > prctile(nostim, 99.7917) | data< prctile(nostim, 0.2083)))];
        n_phases=numel(find(data > prctile(nostim, 99.7917) | data< prctile(nostim, 0.2083)));
        
        if ~isempty(find(data > prctile(nostim, 99.7917) | data< prctile(nostim, 0.2083)))
            phase=find(data > prctile(nostim, 99.7917) | data< prctile(nostim, 0.2083));
        else
            phase=NaN;
        end
        subplot(3,1,rr)
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
        
        %  plot(time,dr,'r.','MarkerSize',5)
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
    end
     p1.OuterPosition= [440,163,201,634];
end


        if type==1
            avg_sup=nanmean(sup,3);
        end
end



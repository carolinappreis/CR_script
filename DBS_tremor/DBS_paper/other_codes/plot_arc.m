function[avg_sup]=plot_arc(out,match_ax,color_b1)
sup=[];
for type=1
    for iii=1:size(out.start_c,1)
        
        feature={'out.mod_amp';'out.mod_freq'}; %%% non stim of freq has not been calculated, uncomment that bit to see frc
        p=type;
        dr=eval([(num2str(feature{type,1})) '{' num2str(iii) ',2}{match_ax(2,iii,1),1}']);
        data=nanmedian(dr);
        
        sup(type,iii)=min(data);
        
        
        nostim=eval(squeeze([num2str(feature{type,1}) '{' num2str(iii) ',1}(match_ax(1,iii,1),:)']));
        n=[];
        n=[n ; data(find(data > prctile(nostim, 99.7917) | data< prctile(nostim, 0.2083)))];
        n_phases=numel(find(data > prctile(nostim, 99.7917) | data< prctile(nostim, 0.2083)));
        
        if ~isempty(find(data > prctile(nostim, 99.7917) | data< prctile(nostim, 0.2083)))
            phase=find(data > prctile(nostim, 99.7917) | data< prctile(nostim, 0.2083));
        else
            phase=NaN;
        end
        
        p1=figure(type+1)
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
%         yline(prctile(nostim,99.7917),'k--','LineWidth',1)
        yline(prctile(nostim,0.2083),'k--','LineWidth',1)
        %  plot(time,dr,'r.','MarkerSize',5)
        if ~isnan(phase)
            plot(time(phase),y2(phase),'.','Color',color_b1(3,:),'MarkerSize',25)

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
        
        %%% all axis
        if type==1
            avg_sup=[nanmean(sup) nanstd(sup)];
            cl=[[0 0 0];[1 0 0];[0 0 1]];
            feature={'out.mod_amp';'out.mod_freq'}; %%% non stim of freq has not been calculated, uncomment that bit to see frc
            p2=figure(type+3)
            subplot(1,4,iii)
            for n=1:3
                dr=eval([(num2str(feature{type,1})) '{' num2str(iii) ',2}{n,1}']);
                data=nanmedian(dr);
                plot(time,data,'LineWidth',2,'Color',cl(n,:))
                hold on 
            end
            dum=nanmedian(out.mod_amp{iii,2}{match_ax(2,iii,1),1});
            xline(time((find(dum==(min(dum))))),'--',{'lock'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','LineWidth',2,'Color',[0.5 0.5 0.5])

            ylim([-1 1])
            ylabel({'Change in tremor severity'})
            xlim([-5 335])
            xticks([0:30:330])
            box('off')
            xlabel({'Stimulation phase (degrees)'})
            set(gca,'FontSize',12)
            set(gca,'FontName','Arial','XTickLabelRotation',45)
            set(p2,'color','w');
            title(sprintf('patient %d',(iii)))
            dumdum= round((mean(mean(out.psi{iii,1}'))),2); dumdum= vpa(dumdum);
            text(180,0.8,sprintf('psi axes = %s',dumdum),'FontSize',10)
             
        end 
        clear dum
    end
    p1.Units = 'centimeters';
    p1.OuterPosition= [10, 10, 35, 10];
    set(p1,'color','w');
    
    p2.Units = 'centimeters';
    p2.OuterPosition= [10, 10, 35, 10];
    set(p2,'color','w');
end

end


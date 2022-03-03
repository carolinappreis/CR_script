function [f1]=modintime_5sec(out,match_ax)
type=1;
for iii=1:size(out.start_c,1)
    
    feature={'out.mod_amp';'out.mod_freq'}; %%% non stim of freq has not been calculated, uncomment that bit to see frc
    dr=eval([(num2str(feature{type,1})) '{' num2str(iii) ',2}{match_ax(2,iii,1),1}']);
    data=nanmedian(dr);
    
    for axi=1 %%% 1=main axis, change to see other axes
        pp1=squeeze(out.seg_zenv{iii,1}{match_ax(2,iii,axi),1}(:,(find(data==(min(data)))),:));
        pp2=squeeze(out.seg_zenv{iii,1}{match_ax(2,iii,axi),1}(:,(find(data==(max(data)))),:));
        dum=find(~isnan(pp1(:,1)));
        for i=1:length(dum)
            pp(1,i,:)=(pp1(i,:)-(mean(pp1(i,1:1000))))./(mean(pp1(i,1:1000)));
            pp(2,i,:)=(pp2(i,:)-(mean(pp2(i,1:1000))))./(mean(pp2(i,1:1000)));
        end
        
        cl=[[1 0 0];[0 0 0];[0 0 1]];
        for modu=1 %%% 1=supression, change to 2 to see amplification
            ppn=squeeze(pp(modu,:,:));
            f1=figure(axi+modu);
            subplot(1,4,iii)
            y2=smooth(nanmean(ppn),500)'; %%% smoothing of half a second applied
            y1=(y2+(smooth(((nanmean(ppn)+nanstd(ppn))./sqrt(size(ppn,1))),500))');
            y3=(y2-(smooth(((nanmean(ppn)+nanstd(ppn))./sqrt(size(ppn,1))),500))');
            time=1:6000;
            patch([time fliplr(time)], [y1 fliplr(y2)],[cl(axi,:)],'FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
            hold on
            patch([time fliplr(time)], [y2 fliplr(y3)],[cl(axi,:)],'FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
            plot(time,y2,'Color',cl(axi,:),'LineWidth',2)
            if modu==1
                ylim([-0.5 0.5])
            else
                ylim([-1 1])
            end
            yticks(-0.5:0.25:0.5)
            xlim([0 6000])
            xticks([0:1000:6000])
            xticklabels({'-1','0','1','2','3','4','5'})
            set(gca,'FontSize',12)
            ylabel('Normalised tremor severity (a.u)')
            xlabel ('time (s)')
            title(sprintf('patient %d',(iii)))
            box('off')
            xline(1000,'--',{'stim ON'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','LineWidth',2,'Color',[0.5 0.5 0.5])
            set(f1,'color','w');
            f1.OuterPosition= [1,100,1000,300];          
        end
    end    
end
end
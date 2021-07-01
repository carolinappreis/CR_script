function[f1]=pt4(out,iii,type,match_ax)

feature={'out.mod_amp';'out.mod_freq'}; %%% non stim of freq has not been calculated, uncomment that bit to see frc
dr=eval([(num2str(feature{type,1})) '{' num2str(iii) ',2}{match_ax(2,iii,1),1}']);
data=nanmedian(dr);

%%%% main axis only
pp1=squeeze(out.seg_zenv{iii,1}{match_ax(2,iii,1),1}(:,(find(data==(max(data)))),:));
dum=find(~isnan(pp1(:,1)));
for i=1:length(dum)
    pp(i,:)=pp1(i,:)./(mean(pp1(i,1:1000)));
end

f1=figure(type+4)
subplot(1,2,1)
y2=smooth(nanmean(pp))';
y1=(y2+smooth(((nanmean(pp)+nanstd(pp))./sqrt(size(pp,1))))');
y3=(y2-smooth(((nanmean(pp)+nanstd(pp))./sqrt(size(pp,1))))');
time=1:6000;
patch([time fliplr(time)], [y1 fliplr(y2)],'r','FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
hold on
patch([time fliplr(time)], [y2 fliplr(y3)],'r','FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
plot(time,y2,'r','LineWidth',2)
xlim([0 6000])
xticks([0:1000:6000])
xticklabels({'-1','0','1','2','3','4','5'})
%     ylim([0 4])
%     yticks(0:1:4)
set(gca,'FontSize',12)
aa=get(gca,'ylim');
ylabel('tremor severity (m/s^2)')
xlabel ('time (s)')
box('off')
xline(1000,'--',{'stim ON'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','LineWidth',2,'Color',[0.5 0.5 0.5]);

%%%  all axes

cl=[[0 0 0];[1 0 0];[0 0 1]];
for n=1:3
    feature={'out.mod_amp';'out.mod_freq'}; %%% non stim of freq has not been calculated, uncomment that bit to see frc
    dr=eval([(num2str(feature{type,1})) '{' num2str(iii) ',2}{match_ax(2,iii,n),1}']);
    data=nanmedian(dr);
    
    pp1=squeeze(out.seg_zenv{iii,1}{match_ax(2,iii,n),1}(:,(find(data==(max(data)))),:));
    dum=find(~isnan(pp1(:,1)));
    for i=1:length(dum)
        pp(i,:)=pp1(i,:)./(mean(pp1(i,1:1000)));
    end
    clear dum
    
    subplot(1,2,2)
    y2=smooth(nanmean(pp))';
    y1=(y2+smooth(((nanmean(pp)+nanstd(pp))./sqrt(size(pp,1))))');
    y3=(y2-smooth(((nanmean(pp)+nanstd(pp))./sqrt(size(pp,1))))');
    time=1:6000;
    patch([time fliplr(time)], [y1 fliplr(y2)],[cl(n,:)],'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off')
    hold on
    patch([time fliplr(time)], [y2 fliplr(y3)],[cl(n,:)],'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off')
    plot(time,y2,'Color',cl(n,:),'LineWidth',1)
    xlim([0 6000])
    xticks([0:1000:6000])
    xticklabels({'-1','0','1','2','3','4','5'})
    %     ylim([0 4])
    %     yticks(0:1:4)
    set(gca,'FontSize',12)
    aa=get(gca,'ylim');
    ylabel('tremor severity (m/s^2)')
    xlabel ('time (s)')
    box('off')
    xline(1000,'--',{'stim ON'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','LineWidth',2,'Color',[0.5 0.5 0.5]) 
end
 set(f1,'color','w');
end

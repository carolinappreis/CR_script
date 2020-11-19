function [stats]=pwelch3(out,stats)

load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','stone')
color_b1={ stone; aegean; blushred};

frange=find(out.F==3):find(out.F==8);
dumm={'out.Pxx_ns_all';'out.Pxx_s_all';'out.Pxx_a_all'};
for i=1:size(dumm,1)
    cr=eval(dumm{i,1});
    pxxrange=cr(:,frange);
    pNSSA(i,:)=max(pxxrange');
    
    
    
    for jo=1:size(cr,1)
        if ~isnan(pxxrange(jo,1))
            dum= frange(find(pxxrange(jo,:) == max(pxxrange(jo,:))));
            fASNS(i,jo)=out.F(dum);
            seg=[dum-3:dum+3];
            sASNS(i,jo,:)=cr(jo,seg); clear seg
        else
            fASNS(i,jo)=NaN;
            sASNS(i,jo,:)=NaN;
        end
    end
    
    clear cr pxxrange dum
end
% kstest(ASNS(1,:)-ASNS(3,:))==1

[p,stats.power_ns_s]=ttest(pNSSA(1,:),pNSSA(2,:))
[p,stats.power_ns_a]=ttest(pNSSA(1,:),pNSSA(3,:))
[p,stats.power_s_a]=ttest(pNSSA(2,:),pNSSA(3,:))
[p,stats.freq_ns_s]=ttest(fASNS(1,:),fASNS(2,:));
[p,stats.freq_ns_a]=ttest(fASNS(1,:),fASNS(3,:));
[p,stats.freq_s_a]=ttest(fASNS(2,:),fASNS(3,:));










close all
f1=figure(18);
subplot(1,3,1)
i=[];
for i=1:3
    dr=squeeze(sASNS(i,:,:));
    y2=median(dr,1);
    y1=prctile(dr,75);
    y3=prctile(dr,25);
    time=1:length(y2);
    patch([time fliplr(time)], [y1 fliplr(y2)],[color_b1{i,1}],'FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
    hold on
    patch([time fliplr(time)], [y2 fliplr(y3)],[color_b1{i,1}],'FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
    plot(time,y2,'LineStyle','-', 'LineWidth',1.5,'Color',color_b1{i,1})
    set(gca,'FontSize',12)
    box('off')
    % ylim([0 1.6])
    % yticks([0:0.4:2])
    xticks([1:1:7])
    xticklabels({'-1.5','-1','-0.5','peak','0.5','1','1.5'})
    clear y2
    hold on
    xlabel('Frequency')
    ylabel('Tremor PSD')
    hold on
    clear y2 y1 y3
end
legend({'NS','SUP','AMP'})
legend('boxoff')
set(gca,'FontSize',12)
legend('boxoff')

subplot(1,3,2)
bar(1:3,[median(fASNS(1,:)) median(fASNS(2,:)) median(fASNS(3,:))],'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.3,'EdgeColor','none')
hold on
subplot(1,3,3)
bar(1:3,[median(pNSSA(1,:)) median(pNSSA(2,:)) median(pNSSA(3,:))],'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.3,'EdgeColor','none')
hold on
for i=1:size(fASNS,2)
    cl= rand(1,3);
    subplot(1,3,2)
    plot(fASNS(1:3,i),'Color',cl,'LineWidth',1)
    hold on
    plot([fASNS(1,i); fASNS(2,i); fASNS(3,i)],'.','MarkerSize',15,'Color',cl)
    % ylim([0 1.6])
    % yticks([0:0.4:2])
    xlim([0 4])
    xticklabels({'NS','SUP','AMP'})
    box('off')
    ylabel('Tremor frequency peak')
    set(gca,'FontSize',12)
    
    subplot(1,3,3)
    plot(pNSSA(1:3,i),'Color',cl,'LineWidth',1)
    hold on
    plot([pNSSA(1,i); pNSSA(2,i); pNSSA(3,i)],'.','MarkerSize',15,'Color',cl)
    ylim([0 1.6])
    yticks([0:0.4:2])
    xlim([0 4])
    xticklabels({'NS','SUP','AMP'})
    box('off')
    ylabel('Tremor power peak')
    set(gca,'FontSize',12)
    
end
f1.Units ='centimeters';
f1.OuterPosition= [5, 5, 30, 20];
set(f1,'color','w');

f2=figure;
subplot(1,2,1)
bar(1:3,[median(fASNS(1,:)) median(fASNS(2,:)) median(fASNS(3,:))],'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.3,'EdgeColor','none')
hold on
subplot(1,2,2)
bar(1:3,[median(pNSSA(1,:)) median(pNSSA(2,:)) median(pNSSA(3,:))],'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.3,'EdgeColor','none')
hold on
for i=1:size(fASNS,2)
    cl= rand(1,3);
    subplot(1,2,1)
    plot(fASNS(1:3,i),'Color',cl,'LineWidth',1)
    hold on
    plot([fASNS(1,i); fASNS(2,i); fASNS(3,i)],'.','MarkerSize',15,'Color',cl)
    % ylim([0 1.6])
    % yticks([0:0.4:2])
    xlim([0 4])
    xticklabels({'NS','SUP','AMP'})
    box('off')
    ylabel('Tremor peak frequency')
    set(gca,'FontSize',12)
    
    subplot(1,2,2)
    plot(pNSSA(1:3,i),'Color',cl,'LineWidth',1)
    hold on
    plot([pNSSA(1,i); pNSSA(2,i); pNSSA(3,i)],'.','MarkerSize',15,'Color',cl)
    % ylim([0 1.6])
    % yticks([0:0.4:2])
    xlim([0 4])
    xticklabels({'NS','SUP','AMP'})
    box('off')
    ylabel('Tremor peak frequency power')
    set(gca,'FontSize',12)
    
end
f2.Units ='centimeters';
% f2.OuterPosition= [257,547,763,331];
set(f2,'color','w');


%% PLS and NS

% % % 
% % % clearvars -except out stats
% % % load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','stone')
% % % color_b1={stone; blushred};
% % % 
% % % frange=find(out.F==3):find(out.F==8);
% % % dumm={'out.Pxx_nsh_all';'out.Pxx_pls_all'};
% % % for i=1:size(dumm,1)
% % %     cr1=eval(dumm{i,1});
% % %     cr=cr1(find(~isnan(cr1(:,1))),:);
% % %     pxxrange=cr(:,frange);
% % %     pASNS(i,:)=max(pxxrange');
% % %     
% % %     for jo=1:size(cr,1)
% % %         if ~isnan(pxxrange(jo,1))
% % %             dum= frange(find(pxxrange(jo,:) == max(pxxrange(jo,:))));
% % %             fASNS(i,jo)=out.F(dum);
% % %             seg=[dum-3:dum+3];
% % %             sASNS(i,jo,:)=cr(jo,seg); clear seg
% % %         else
% % %             fASNS(i,jo)=NaN;
% % %             sASNS(i,jo,:)=NaN;
% % %         end
% % %     end
% % %     
% % %     clear cr1 cr pxxrange dum
% % % end
% % % 
% % % [p,stats.power_ns_pls]=ttest(pASNS(1,:),pASNS(2,:))
% % % [p,stats.freq_ns_pls]=ttest(fASNS(1,:),fASNS(2,:))
% % % 
% % % 
% % % f1=figure(19);
% % % subplot(1,3,1)
% % % for i=1:2
% % %     dr=squeeze(sASNS(i,:,:));
% % %     y2=nanmedian(dr,1);
% % %     y1=prctile(dr,75);
% % %     y3=prctile(dr,25);
% % %     time=1:length(y2);
% % %     patch([time fliplr(time)], [y1 fliplr(y2)],[color_b1{i,1}],'FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
% % %     hold on
% % %     patch([time fliplr(time)], [y2 fliplr(y3)],[color_b1{i,1}],'FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
% % %     plot(time,y2,'LineStyle','-', 'LineWidth',1.5,'Color',color_b1{i,1})
% % %     set(gca,'FontSize',12)
% % %     box('off')
% % %     % ylim([0 1.6])
% % %     % yticks([0:0.4:2])
% % %     xticks([1:1:7])
% % %     xticklabels({'-1.5','-1','-0.5','peak','0.5','1','1.5'})
% % %     clear y2
% % %     hold on
% % %     xlabel('Frequency')
% % %     ylabel('Tremor PSD')
% % %     hold on
% % %     clear y2 y1 y3 dr
% % % end
% % % legend({'NS','PLS'})
% % % legend('boxoff')
% % % set(gca,'FontSize',12)
% % % legend('boxoff')
% % % 
% % % subplot(1,3,2)
% % % bar(1:2,[median(fASNS(1,:))  median(fASNS(2,:))],'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.3,'EdgeColor','none')
% % % hold on
% % % subplot(1,3,3)
% % % bar(1:2,[ median(pASNS(1,:)) median(pASNS(2,:))],'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.3,'EdgeColor','none')
% % % hold on
% % % for i=1:size(fASNS,2)
% % %     cl= rand(1,3);
% % %     subplot(1,3,2)
% % %     plot(fASNS(1:2,i),'Color',cl,'LineWidth',1)
% % %     hold on
% % %     plot([ fASNS(1,i); fASNS(2,i)],'.','MarkerSize',15,'Color',cl)
% % %     % ylim([0 1.6])
% % %     % yticks([0:0.4:2])
% % %     xlim([0 3])
% % %     xticklabels({'NS','PLS'})
% % %     box('off')
% % %     ylabel('Tremor frequency peak')
% % %     set(gca,'FontSize',12)
% % %     
% % %     subplot(1,3,3)
% % %     plot(pASNS(1:2,i),'Color',cl,'LineWidth',1)
% % %     hold on
% % %     plot([ pASNS(1,i); pASNS(2,i)],'.','MarkerSize',15,'Color',cl)
% % %     % ylim([0 1.6])
% % %     % yticks([0:0.4:2])
% % %     xlim([0 3])
% % %     xticklabels({'NS','PLS'})
% % %     box('off')
% % %     ylabel('Tremor power peak')
% % %     set(gca,'FontSize',12)
% % %     
% % % end
% % % f1.Units ='centimeters';
% % % f1.OuterPosition= [5, 5, 30, 20];
% % % set(f1,'color','w');

end






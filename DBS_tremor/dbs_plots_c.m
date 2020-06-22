function [r,p]=plots_c(out)
m_ax=1;% change if the main axis is not always 1 - replace by array of main axes
load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','stone','squash','sapphire','azure');
color_b1=[blushred; aegean; stone];


%%  plot arc's main axis with significant threshold\
for iii=1:10
    %%% for type=1:2
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
    
    figure(7+type)
    subplot(2,5,iii)
    y2=data;
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
%     yline(prctile(nostim,99.7917),'k--','LineWidth',1)
%     yline(prctile(nostim,0.2083),'k--','LineWidth',1)
    ylim([-1 1])
    xlim([-5 335])
    xticks([0:30:330])
    box('off')
    ylabel({'Change in tremor severity'})
    xlabel({'Stimulation phase (degrees)'})
    set(gca,'FontSize',9)
    set(gca,'FontName','Arial','XTickLabelRotation',45)
end


%% plot median split-arc main axis no thresholds

m_ax=1;% change if the main axis is not always 1 - replace by array of main axes
load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','stone','squash','sapphire','azure');
color_b1=[blushred; aegean; stone];

for iii=1:10
    f1=figure(10)
    subplot(2,5,iii)
    pp=[sapphire;azure];
    
    y=[];clear i
    for i=1:2
        datas(i,:)=eval(['nanmedian(out.arc' num2str(i) '{' num2str(iii) ',2}{m_ax,1})']);
        raw(i,:)=repmat(datas(i,:),1,3);
        for ti=size(datas,2)+1:size(datas,2)*2
            smo_m(i,ti-12)=sum(raw(i,(ti-1:ti+1)))./length(raw(i,(ti-1:ti+1)));
        end
        nsdatas(i,:)=squeeze(eval(['out.ns_ms{' num2str(iii) ',1}(' num2str(i) ',:)']))';
        if ~isempty(find(datas(i,:) > prctile(nsdatas(i,:), 99.7917) | datas(i,:)< prctile(nsdatas(i,:), 0.2083)))
            pha{i,1}=find(datas(i,:) > prctile(nsdatas(i,:), 99.7917) | datas(i,:)< prctile(nsdatas(i,:), 0.2083));
        else
            pha{i,1}=NaN;
        end
        dum=datas;
        %         dum=smo_m;
        %         yline(0,'LineWidth',0.5,'Color', [0.5 0.5 0.5])
        hold on
        stem(1:12,dum(i,:),'.', 'LineWidth',2,'MarkerSize',10,'Color',pp(i,:))
        modms(iii,i)=mean(dum(i,:));
        %                 yline(prctile(nsdatas(i,:),99.7917),'--','Color',pp(i,:))
        %                 yline(prctile(nsdatas(i,:),0.2083),'--','Color',pp(i,:))
        %         ll=round(max([abs(min(min(dum))) abs(max(max(dum)))]),1)+0.1;
        %         ylim([-ll ll]);
        xticks([1:12])
        xticklabels([0:30:330])
    end
    
    for i=1:2
        plot(1:12,dum(i,:),'.', 'LineWidth',1,'MarkerSize',10,'Color',pp(i,:))
        if ~isnan(pha{i,1})
            plot(pha{i,1},dum(i,pha{i,1}),'.','MarkerSize',15,'Color',blushred)
        end
    end
    
    
    box('off')
    ylabel({'Change in tremor severity'})
    xlabel({'Stimulation phase (degrees)'})
    set(gca,'FontSize',12)
    set(gca,'FontName','Arial','XTickLabelRotation',45)
     set(f1,'color','w');
    
    clear i
end

    f1=figure(31)
    b=bar(median(abs(modms)),'FaceColor','flat','FaceAlpha',.6,'EdgeColor','none','BarWidth',1)
    b.CData(1,:) = pp(1,:);
    b.CData(2,:) = pp(2,:);
    hold on
    plot(abs(modms)','k.')
    xticklabels({'low amp','high amp'})
    xlabel({'tremor state'})
    ylabel({'mean magnitude';'of modulation'})
    box('off')
    ylim([0 0.8 ])
    set(f1,'color','w');
    set(gca,'FontSize',12);
     ttest(abs(modms(:,1)),abs(modms(:,2)))
%     ttest(modms(:,1),modms(:,2))
    
%     er = errorbar(1:2,median(modms),prctile(modms,5),prctile(modms,95));
%     er.Color = [0 0 0];
%     er.LineStyle = 'none';

%% arcs/arc_smothed all axes

m_ax=1;% change if the main axis is not always 1 - replace by array of main axes
load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','stone','squash','sapphire','azure','ax2','ax3');
color_b1=[stone;ax2;ax3];


for iii=1:10
    t_sig=[];
    degs=[0:30:330];
    for ax=1:3
        raw_s=squeeze(median(out.change_c{iii,2}{ax,1}));
        enve=squeeze(median(out.end_env{iii,2}{ax,1}));
        raw_ns=squeeze(out.change_c{iii,1}(ax,:));
        sig=find(raw_s > prctile(raw_ns, 99.7917) | raw_s< prctile(raw_ns, 0.2083));
        t_sig=[t_sig degs(sig)];
        
        raw_s1=repmat(raw_s,1,3);
        for i=size(raw_s,2)+1:size(raw_s,2)*2
            smo_s(1,i-12)=sum(raw_s1(1,(i-1:i+1)))./length(raw_s1(1,(i-1:i+1)));
        end
        
        %choice=smo_s;
        choice=raw_s;
        a_arc(1,ax)=abs(sum(raw_s));
        a_env(1,ax)=abs(sum(enve));
        a_abs(1,ax)=mean(enve);
        f1=figure(11);
        subplot(2,5,iii)
        plot(1:12,choice,'Color',color_b1(ax,:),'LineWidth',1.5)
        hold on
        %         if ~isempty(sig)
        %             plot(sig,choice(1,sig),'*','MarkerSize',5,'Color','r')
        %         end
        xticks([1:12])
        xticklabels([0:30:330])
        %         ylim([-1.5 1.5])
        ylim([-1.5 1.5])
        ylabel('Change in tremor severity')
        xlabel('Stimulation phase (degrees)')
        set(gca,'XTickLabelRotation',45)
        box('off')
        set(f1,'color','w');
        clear sig raw_s raw_ns smo_s
    end
    tall_sig{iii,1}=t_sig;
    ax_mod(iii,:)=a_arc./sum(a_arc)*100;
    %     amp_mod(iii,:)=a_env./sum(a_env)*100;
    amp_abs(iii,:)=a_abs;
    %     ax_abs(iii,:)=a_arc;
    ax_max(iii,1)=find(ax_mod(iii,:)==(max(ax_mod(iii,:))));
    m_max(iii,1)=find(amp_abs(iii,:)==(max(amp_abs(iii,:))));
    
    t_s=tall_sig(~(cellfun('isempty',tall_sig)));
    dum=0.15707963267949;
    f1=figure(12) %%% error - not plotting first subplot correctly
    for k=1:size(t_s,1)
        subplot(1,size(t_s,1),k)
        polarhistogram(deg2rad(t_s{k,1}),30,'FaceColor','red','FaceAlpha',.8,'EdgeColor','none','BinWidth',dum)
        ax=gca;
        ax.RLim=[0 3];
        ax.RTick=[0 1 2 3];
        ax.GridLineStyle = '--';
        ax.RGrid='off'
        ax.GridAlpha = 0.3;
        ax.Box = 'off';
        clear x
    end
    k=1
    subplot(1,size(t_s,1),k)
    polarhistogram(deg2rad(t_s{k,1}-pi),30,'FaceColor','red','FaceAlpha',.8,'EdgeColor','none','BinWidth',dum)
    ax=gca;
    ax.RLim=[0 3];
    ax.RTick=[0 1 2 3];
    ax.GridLineStyle = '--';
    ax.RGrid='off'
    ax.GridAlpha = 0.3;
    ax.Box = 'off';
    clear x
    set(f1,'color','w');
end

% figure(13)
% subplot(2,2,1)
% bar(ax_mod)
% xlabel('patients')
% ylabel('tremor change(%)')
% legend('main axis','axis2','axis3')
% ylim([0 100])
% legend('boxoff')
% box('off')
% subplot(2,2,2)
% bar([numel(find(ax_max==1)) numel(find(ax_max==2|ax_max==3))])
% xticklabels({'main axis','other axes'})
% set(gca,'FontName','Arial','XTickLabelRotation',45)
% ylabel('patients')
% ylim([0 10])
% box('off')
% subplot(2,2,3)
% bar(amp_mod)
% xlabel('patients')
% ylabel('tremor severity(%)')
% legend('main axis','axis2','axis3')
% ylim([0 100])
% legend('boxoff')
% box('off')
% subplot(2,2,4)
% bar([numel(find(m_max==1)) numel(find(m_max==2|m_max==3))])
% ylim([0 10])
% xticklabels({'main axis','other axes'})
% set(gca,'FontName','Arial','XTickLabelRotation',45)
% ylabel('patients')
% box('off')

f1=figure(13)
subplot(2,2,1)
b=bar(ax_mod,'FaceAlpha',[0.75],'EdgeColor','none');
b(1).FaceColor=stone;
b(2).FaceColor=ax2;
b(3).FaceColor=ax3;
xlabel('patients')
ylabel({'change in tremor' ; 'severity (%)'})
legend('main axis','axis2','axis3')
legend('boxoff')
box('off')
set(gca,'FontSize',12);
subplot(2,2,2)
b1=bar([numel(find(ax_max==1)) numel(find(ax_max==2|ax_max==3))],'EdgeColor','none','FaceColor','flat');
b1.CData(1,:) = stone;
b1.CData(2,:) = [1.00,0.54,0.16];
xticklabels({'main axis','other axes'})
% set(gca,'FontName','Arial','XTickLabelRotation',45)
ylabel('patients')
ylim([0 10])
box('off')
set(gca,'FontSize',12);
subplot(2,2,3)
b=bar(amp_abs,'FaceAlpha',[0.75],'EdgeColor','none')
b(1).FaceColor=stone;
b(2).FaceColor=ax2;
b(3).FaceColor=ax3;
xlabel('patients')
ylabel({'tremor' ;' severity (m/s^2)'})
legend('main axis','axis2','axis3')
legend('boxoff')
box('off')
set(gca,'FontSize',12);
subplot(2,2,4)
b1=bar([numel(find(m_max==1)) numel(find(m_max==2|m_max==3))],'EdgeColor','none','FaceColor','flat')
b1.CData(1,:) = stone;
b1.CData(2,:) = [1.00,0.54,0.16];
ylim([0 10])
xticklabels({'main axis','other axes'})
% set(gca,'FontName','Arial','XTickLabelRotation',45)
ylabel('patients')
box('off')
set(gca,'FontSize',12);
set(f1,'color','w');


%% gaussian fit to tremor properties & correlation with arc features
load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','stone');
cl=stone;
m_ax=1;
for iii=1:10
    
    f1=figure(14)
    subplot(2,5,iii)
    y=out.amp_n_bins(iii,:);
    width(iii,:)=sum(~isnan(y));
    x=out.bins(1:end-2);
    
    bar(x,y,'FaceColor',cl,'EdgeColor',cl)
    hold on
    [rsg,rsg_g,rsg_o]=gauss_fit2(x,y)
    ylim([0 1.5])
    xticks([1:2:14])
    %     xticklabels({'2','3','4','5','6','7','8'});
    ylabel({'absolute amplitude \uV^2'})
    xlabel('frequency(Hz)')
    %     set(gca,'FontSize',12)
    box('off')
    legend('off')
    cv(iii,:)=rsg.c./rsg.b;
    clear y rsg rsg_o rsg_g
    
    f1.Units = 'centimeters';
    f1.OuterPosition= [10, 10, 25, 10];
    set(f1,'color','w');
    
    
    d1=nanmedian(squeeze(out.change_c{iii,2}{1,m_ax}));
    arc.value(iii,1,1)=nanmean(abs(d1))*100;
    arc.value(iii,2,1)=(max(d1)-min(d1));
    if min(d1)<0
        arc.value(iii,3,1)=-(min(d1)*100);
    else
        arc.value(iii,3,1)=NaN;
    end
    %     arc.value(iii,4,1)=max(d1)*100;
    arc.label={'mean effect';'effect range';'supressive effect'};
    clear d1
end

for ii=1:size(arc.value,2)
    
    f1=figure(15)
    subplot(1,size(arc.value,2),ii)
    
    dum=find(~isnan((arc.value(:,ii))));
    x=cv(dum);
    %     x=width(dum);
    y=arc.value(dum,ii);
    
    [fitresult] = fit( x, y, 'poly1' );
    h=plot(fitresult,x,y);
    
    hold on
    plot(x,y,'k.','MarkerSize',10,'HandleVisibility','off');
    legend('off')
    [r(1,ii),p(1,ii)]=corrcoef(x,y);
    
    %     legend(h,['r=',num2str(r)],'box','off')
    box('off')
    
    xlabel('coefficient of variation')
    ylabel(sprintf( arc.label{ii,1}))
    
    
    set(gca,'FontSize',12)
    f1.Units ='centimeters';
    f1.OuterPosition= [10, 10, 8, 10];
    set(f1,'color','w');
    
    %     y2=plot(x(sig),y(sig),'ko','MarkerSize',8);
    
end

%% PLS vs NS

% % figure(16)
% % yu=find(~isnan(out.pls3(:,1)));
% % for o=1:length(yu)
% %     subplot(1,5,o)
% %     plot(out.ns3(yu(o),:),'.-','Color',stone,'LineWidth',2,'MarkerSize',10)
% %     hold on
% %     plot(out.pls3(yu(o),:),'.-','Color',blushred,'LineWidth',2,'MarkerSize',10)
% %     xlim([0.5 3.5])
% %     box('off')
% %     ylabel({'Tremor severity';'(zscore)'})
% %     xticklabels({'start','middle','end'})
% %     xtickangle(45)
% % end
% %
% % figure(17)
% % plot(out.ns3(yu,:)','o','Color',stone,'MarkerSize',5)
% % hold on
% % plot(out.pls3(yu,:)','o','Color',blushred,'MarkerSize',5)
% % plot(median(out.ns3(yu,:)),'.-','Color',stone,'LineWidth',2,'MarkerSize',15)
% % plot(median(out.pls3(yu,:)),'.-','Color',blushred,'LineWidth',2,'MarkerSize',15)
% % xlim([0.5 3.5])
% % xticks([1 2 3])
% % box('off')
% % ylabel({'Tremor severity';'(zscore)'})
% % xticklabels({'start','middle','end'})
% % xtickangle(45)
% % set(gca,'FontSize',12)
% %
% % [f,h]=ttest(out.ns3(yu,:),out.pls3(yu,:))
end


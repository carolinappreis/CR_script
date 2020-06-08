function [r,p]=plots_c(out)
for iii=1:10
    m_ax=1;% change if the main axis is not always 1 - replace by array of main axes
    
    load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','stone','squash','sapphire','azure');
    color_b1=[blushred; aegean; stone];
    % load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/clust_master_output.mat')
    
    
    %%%%  plot arc's main axis with significant threshold
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
    
    
    
    %%%% plot ampsplit-arc main axis no thresholds
    
    load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/NS_all.mat','ns_ba_split');
    
    figure(10)
    subplot(2,5,iii)
    pp=[azure;sapphire];
    
    y=[];clear i
    for i=1:2
        datas(i,:)=eval(['nanmedian(out.arc' num2str(i) '{' num2str(iii) ',2}{m_ax,1})']);
        raw(i,:)=repmat(datas(i,:),1,3);
        for ti=size(datas,2)+1:size(datas,2)*2
            smo_m(i,ti-12)=sum(raw(i,(ti-1:ti+1)))./length(raw(i,(ti-1:ti+1)));
        end
        nsdatas(i,:)=squeeze(eval(['ns_ba_split(' num2str(iii) ',' num2str(i) ',:)']))';
        if ~isempty(find(datas(i,:) > prctile(nsdatas(i,:), 99.7917) | datas(i,:)< prctile(nsdatas(i,:), 0.2083)))
            pha{i,1}=find(datas(i,:) > prctile(nsdatas(i,:), 99.7917) | datas(i,:)< prctile(nsdatas(i,:), 0.2083));
        else
            pha{i,1}=NaN;
        end
        %         dum=datas;
        dum=smo_m;
        %         yline(0,'LineWidth',0.5,'Color', [0.5 0.5 0.5])
        hold on
        stem(1:12,dum(i,:),'.', 'LineWidth',2,'MarkerSize',10,'Color',pp(i,:))
        
        %         yline(prctile(nsdatas(i,:),99.7917),'--','Color',pp(i,:))
        %         yline(prctile(nsdatas(i,:),0.2083),'--','Color',pp(i,:))
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
    
    clear i
    
    %%% arc_smothed all axes
    for ax=1:3
        raw_s=squeeze(median(out.change_c{iii,2}{ax,1}));
        raw_ns=squeeze(out.change_c{iii,1}(ax,:));
        sig=find(raw_s > prctile(raw_ns, 99.7917) | raw_s< prctile(raw_ns, 0.2083));
        
        raw_s1=repmat(raw_s,1,3);
        for i=size(raw_s,2)+1:size(raw_s,2)*2
            smo_s(1,i-12)=sum(raw_s1(1,(i-1:i+1)))./length(raw_s1(1,(i-1:i+1)));
        end
        
        choice=smo_s;
        
        figure(11);
        subplot(2,5,iii)
        plot(1:12,choice,'Color',color_b1(ax,:),'LineWidth',1.5)
        hold on
        if ~isempty(sig)
            plot(sig,choice(1,sig),'*','MarkerSize',5,'Color','r')
        end
        xticks([1:12])
        xticklabels([0:30:330])
        %         ylim([-1.5 1.5])
        ylim([-1 1])
        ylabel('Change in tremor severity')
        xlabel('Stimulation phase (degrees)')
        set(gca,'XTickLabelRotation',45)
        box('off')
        clear sig raw_s raw_ns smo_s
    end
    
    %%%% gaussian fit to tremor properties & correlation with arc features
    
    load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','stone');
    cl=stone;
    
    
    f1=figure(12)
    subplot(2,5,iii)
    y=out.amp_n_bins(iii,:);
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
    %%%%------
    
end


for ii=1:size(arc.value,2)
    
    f1=figure(13)
    subplot(1,size(arc.value,2),ii)
    
    dum=find(~isnan((arc.value(:,ii))));
    x=cv(dum);
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

%%% PLS vs NS
yu=find(~isnan(out.pls4(:,1)));
for o=1:length(yu)
    subplot(1,5,o)
    plot(out.ns4(yu(o),:),'k.-','LineWidth',2)
    hold on
    plot(out.pls4(yu(o),:),'r.-','LineWidth',2)
    xlim([0.5 4.5])
    box('off')
    ylabel({'tremor severity (zscore)'})
    %%% add xticklabels
end



end
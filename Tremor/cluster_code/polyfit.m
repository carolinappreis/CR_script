function[out]=polyfit(out)
for iii=1:size(out.start_c,1)
    m_ax=1;
    f1=figure(7)
    if size(out.start_c,1)==4
    subplot(1,4,iii)
    else
    subplot(2,5,iii)
    end
    x = 1:12;
    y = squeeze(nanmedian(out.change_c{iii,2}{1,m_ax}));
    bar(x,y,'FaceColor',[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.5)
    box('off')
    hold on
    y=y';x=x';
    xData=x(find(~isnan(y))); yData=y(~isnan(y));
    mdp1=fitlm(x,y,'poly1')
    mdp2=fitlm(x,y,'poly2')
    mdp3=fitlm(x,y,'poly3')
    mdp4=fitlm(x,y,'poly4')
    mdk= fitlm(x,y,'constant');
    plot(mdp1.Fitted,'y','LineWidth',1.5)
    plot(mdp2.Fitted,'b','LineWidth',1.5)
    plot(mdp3.Fitted,'r','LineWidth',1.5)
    plot(mdp4.Fitted,'g','LineWidth',1.5)
    plot(mdk.Fitted,'k','LineWidth',1.5)
    xticks([1:12])
    xticklabels([0:30:330])    
    box('off')
    ylabel({'Change in tremor severity'},'FontSize',12)
    xlabel({'Stimulation phase (degrees)'},'FontSize',12)
    set(gca,'FontName','Arial','XTickLabelRotation',50)
    set(f1,'color','w');
    title(sprintf('patient %d',(iii)))

    
    rsqr(iii,:)=[ mdk.Rsquared.Adjusted mdp1.Rsquared.Adjusted mdp2.Rsquared.Adjusted mdp3.Rsquared.Adjusted  mdp4.Rsquared.Adjusted ];
    
    F_stat=[ coefTest(mdk) coefTest(mdp1) coefTest(mdp2) coefTest(mdp3)  coefTest(mdp4)];
    [r_val r_idx]=sort(rsqr(iii,:),'descend');
    out.win_m(iii,:)=[r_idx(1) r_val(1) F_stat(r_idx(1))];
%     
%     F_stat=[ coefTest(mdp1) coefTest(mdp2) coefTest(mdp3)  coefTest(mdp4)];
%     [r_val r_idx]=sort(F_stat,'ascend');
%     out.win_m(iii,:)=[r_idx(1) F_stat(r_idx(1))];
    
    
    %     mcomp (iii,:)=[mdk.ModelCriterion.AICc mdp1.ModelCriterion.AICc mdp2.ModelCriterion.AICc mdp3.ModelCriterion.AICc mdp4.ModelCriterion.AICc];
    %     [s_val s_idx]=sort(mcomp(iii,:),'ascend');
    %     win(iii,1)=s_idx(1);
    %     dif_aic(iii,:)=[(abs(s_val(1)-s_val(2)))>2 (abs(s_val(1)-s_val(3)))>2 (abs(s_val(1)-s_val(4)))>2 (abs(s_val(1)-s_val(5)))>2];
end
 f1.OuterPosition= [1,339,1440,539];

end
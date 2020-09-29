function [win]=aid_fx(iii,uu,x,y)
win=[];
xData=x(find(~isnan(y))); yData=y(~isnan(y));
mdp1=fitlm(x,y,'poly1')
mdp2=fitlm(x,y,'poly2')
mdp3=fitlm(x,y,'poly3')
mdp4=fitlm(x,y,'poly4')
mdk= fitlm(x,y,'constant');
% %     bar(x,y,'FaceColor',[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.5)
% %     hold on
% %     box('off')
% %         plot(mdp1.Fitted,'y','LineWidth',1.5)
% %         plot(mdp2.Fitted,'b','LineWidth',1.5)
% %         plot(mdp3.Fitted,'r','LineWidth',1.5)
% %         plot(mdp4.Fitted,'g','LineWidth',1.5)
% %         plot(mdk.Fitted,'k','LineWidth',1.5)
%         xticks([1:12])
%         xticklabels([0:30:330])
%         box('off')
%         ylabel({'Change in tremor severity'},'FontSize',12)
%         xlabel({'Stimulation phase (degrees)'},'FontSize',12)
%         set(gca,'FontName','Arial','XTickLabelRotation',50)
%         set(f1,'color','w');


rsqr(uu,:)=[ mdk.Rsquared.Adjusted mdp1.Rsquared.Adjusted mdp2.Rsquared.Adjusted mdp3.Rsquared.Adjusted  mdp4.Rsquared.Adjusted ];

%   F_stat=[ coefTest(mdk) coefTest(mdp1) coefTest(mdp2) coefTest(mdp3)  coefTest(mdp4)];
[r_val r_idx]=sort(rsqr(uu,:),'descend');
win=[r_idx(1) r_val(1)];

end
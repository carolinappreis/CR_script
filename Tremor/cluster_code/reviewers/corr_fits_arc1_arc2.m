% clear; close
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cluster_out_mc.mat','out');
close; clearvars -except out
for iii=1:size(out.arc1,1)
    for t=1:2
        feature={'out.arc1';'out.arc2'};
        dr=eval([(num2str(feature{t,1})) '{' num2str(iii) ',2}{1,1}']);
        data(iii,t,:)=nanmedian(dr);
        x=1:size(data,3);
        y=smooth(squeeze(data(iii,t,:)));
        hold on
        mdp1=fitlm(x,y,'poly6');
        data_fit(iii,t,:)=mdp1.Fitted';
        
        fig=figure(1);
        cl=[[0 0 0];[1 0 0]];
        subplot(2,size(out.arc1,1)/2,iii)
        plot(x,squeeze(data(iii,t,:)),'Color',cl(t,:))
        hold on
        plot(x,squeeze(data_fit(iii,t,:)),'-','Color',cl(t,:),'LineWidth',1.5)
        xlim([0 12])
        xticks(1:1:12)
        ylim([-1 1])
        box('off')
        clear mdp1 dr x y  r p
    end
    
%    [r,p]=corrcoef(squeeze(data_fit(iii,1,:)),squeeze(data_fit(iii,2,:)));
    [r,p]=corrcoef(squeeze(data(iii,1,:)),squeeze(data(iii,2,:)));

    cor(iii,:)=[r p];
    if cor(iii,2)<0.05
        filename=['corr=',num2str(round(cor(iii,1),2)),'*'];
    else
        filename=['corr=',num2str(round(cor(iii,1),2))];
    end
    text(10,0.9,filename);
    clearvars -except out iii data data_fit fig cir filename
    
end

set(fig,'color','w');
fig.OuterPosition= [1,339,1440,539];

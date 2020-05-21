
clear all; close all;
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/new_cluster_CR')
load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','squash');
cl=blushred;
tt1=NaN(10,3,12);
for i=1:10
    for ii=1:3
        tt1(i,ii,:)=nanmedian(tt1_all{i,ii})
    end
end


n=[]; main_clust=[1 1 1 1 1 1 1 1 1 1];
for i=1:size(tt1,1)
%     f1 = figure(i)
    for axis = 1:size(tt1, 2)
        
%         data=squeeze(tt1(i, axis,:));
%         n=[n ; data(find(data > prctile(nostim(i, axis, :), 99.7917) | data< prctile(nostim(i, axis, :), 0.2083)))];
%         n_phases(i,axis)=numel(find(data > prctile(nostim(i, axis, :), 99.7917) | data< prctile(nostim(i, axis, :), 0.2083)));
%         if ~isempty(find(data > prctile(nostim(i, axis, :), 99.7917) | data< prctile(nostim(i, axis, :), 0.2083)))
%             phase{i,axis}=find(data > prctile(nostim(i, axis, :), 99.7917) | data< prctile(nostim(i, axis, :), 0.2083));
%         else
%             phase{i,axis}=NaN;
%         end
        
        %         subplot(1, 3, axis)
        %         bar(0:30:330, squeeze(tt1(i, axis,:)),'FaceColor',cl,'EdgeColor',cl)
        %         hold on
        %         yline(prctile(nostim(i, axis, :), 99.7917),'k--','LineWidth',1)
        %         yline(prctile(nostim(i, axis, :), 0.2083),'k--','LineWidth',1)
        %         %                 ylim([-1 1])
        %         ylabel('Change in tremor severity')
        %         xlabel('Stimulation phase (degrees)')
        %         set(gca,'XTickLabelRotation',45)
        %         box('off')
        %         f1.Units = 'centimeters';
        %         f1.OuterPosition= [10, 10, 30, 8];
        %         set(f1,'color','w');
        %         filename=['NEW_clust',num2str(i),'.png'];
        % %         saveas(gcf,filename)
        
       f1 = figure(1)
        subplot(2,5,i)
        axis=1
        bar(0:30:330, squeeze(tt1(i, axis,:)),'FaceColor',cl,'EdgeColor',cl)
        hold on
        yline(prctile(nostim(i, axis, :), 99.7917),'k--','LineWidth',1)
        yline(prctile(nostim(i, axis, :), 0.2083),'k--','LineWidth',1)
        ylim([-0.8 0.8])
        box('off')
        ylabel({'Change in tremor severity'},'FontSize',12)
        xlabel({'Stimulation phase (degrees)'},'FontSize',12)
        set(gca,'FontName','Arial','XTickLabelRotation',50)
%         f2=figure(2)
%         subplot(2,5,i)
%         histogram(nostim(i,axis,:))


        f1=figure(2);
        subplot(2,5,i)
        x = 1:12;
        dr=squeeze(tt1_all{i,1});
        y=nanmedian(dr,1);
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
        plot(mdp1.Fitted,'k','LineWidth',1.5)
        plot(mdp2.Fitted,'b','LineWidth',1.5)
        plot(mdp3.Fitted,'r','LineWidth',1.5)
        plot(mdp4.Fitted,'g','LineWidth',1.5)
        %plot(mdk.Fitted,'k','LineWidth',1.5)
        %     xticks([1:12])
        %     xticklabels([0:30:330])
        %     a = get(gca,'XTickLabel');
        %     set(gca,'XTickLabel',a,'FontSize',5)
        %     b = get(gca,'YTickLabel');
        %     set(gca,'YTickLabel',b,'FontSize',10)
        %
        %     box('off')
        %     ylabel({'Change in tremor severity'},'FontSize',12)
        %     xlabel({'Stimulation phase (degrees)'},'FontSize',12)
        %     set(gca,'FontName','Arial','XTickLabelRotation',50)
        %
    
        rsqr(i,:)=[ mdp1.Rsquared.Adjusted mdp2.Rsquared.Adjusted mdp3.Rsquared.Adjusted  mdp4.Rsquared.Adjusted ];
        F_stat(i,:)=[ coefTest(mdp1) coefTest(mdp2)  coefTest(mdp3)  coefTest(mdp4)];
        [r_val r_idx]=sort(rsqr(i,:),'descend');
        win_m(i,:)=[r_idx(1) r_val(1) F_stat(i,r_idx(1))];
    
    
        mcomp (i,:)=[mdk.ModelCriterion.AICc mdp1.ModelCriterion.AICc mdp2.ModelCriterion.AICc mdp3.ModelCriterion.AICc mdp4.ModelCriterion.AICc];
        [s_val s_idx]=sort(mcomp(i,:),'ascend');
        win(i,1)=s_idx(1);
        dif_aic(i,:)=[(abs(s_val(1)-s_val(2)))>2 (abs(s_val(1)-s_val(3)))>2 (abs(s_val(1)-s_val(4)))>2 (abs(s_val(1)-s_val(5)))>2];
%         
        
    end
end






% sig_calc=size(tt1,1)*size(tt1,2)*size(tt1,3)*0.05;
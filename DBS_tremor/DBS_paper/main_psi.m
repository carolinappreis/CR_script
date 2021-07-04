function[p1]=main_psi(out,color_b1)

for type=1:2
    p1=figure(type);
    for iii=1:size(out.start_c,1)
        feature={'out.axe_sum';'out.psi'}; %%% non stim of freq has not been calculated, uncomment that bit to see frc
        dr=eval((num2str(feature{type,1})));
        if type==1
            subplot(1,4,iii)
            bar(dr(iii,:),'EdgeColor','none','FaceColor',color_b1(1,:),'FaceAlpha',0.5);
            ylim([0 80])
            yticks(0:20:80)
            xticklabels({'Z','Y','X'})
            ylabel('trials')
            box('off')
            set(gca,'FontSize',14)
            
        elseif type==2
            subplot(1,4,iii)
            err=[std(dr{iii,1}(1,:)) std(dr{iii,1}(2,:))];
            x=[mean(dr{iii,1}(1,:)) mean(dr{iii,1}(2,:))];
            bar([mean(dr{iii,1}(1,:)) mean(dr{iii,1}(2,:))],'EdgeColor','none','FaceColor',color_b1(1,:),'FaceAlpha',0.5)
            hold on
            errorbar([1:2],x,err,'.','Color',color_b1(1,:),'LineWidth',2)
            ylim([0 1.2])
            yticks(0:0.5:1)
            xticklabels({['z,y'],['z,x']}) 
            xlabel ('Axes')
            ylabel('PSI')
            box('off')
            set(gca,'FontSize',14)
            
        end
    end
    
    p1.OuterPosition= [1,100,1000,300];
    set(p1,'color','w');
    
end
end



% % label={'Normalised Power' ; 'Normalised Power';'Mag-Squared Coherence';};
% % 
% %  ylabel(sprintf(label{ii,1}));

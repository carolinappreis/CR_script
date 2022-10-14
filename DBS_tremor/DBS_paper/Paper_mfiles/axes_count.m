function[p1]=axes_count(out)

type=1;
    p1=figure(type);
    for iii=1:size(out.start_c,1)
        feature={'out.axe_sum';'out.psi'}; %%% non stim of freq has not been calculated, uncomment that bit to see frc
        dr=eval((num2str(feature{type,1})));
        total=sum(dr(iii,:));
        share=[(dr(iii,1)*100)/total (dr(iii,2)*100)/total (dr(iii,3)*100)/total];
            subplot(1,4,iii)
            pie(share)
            box('off')
            title(sprintf('patient %d',(iii)))

            set(gca,'FontSize',14)
            clear total share
    end

set(p1,'color','w')
legend('Z','Y','X')
legend('boxoff')
end



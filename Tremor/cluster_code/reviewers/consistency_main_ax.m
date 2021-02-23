clear ; close 
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cluster_out_mc.mat');

for sub=1:10
    for pp=1:12
        for ii=1:3
            ax(:,ii)=out.end_env{sub,2}{ii,1}(:,pp);
        end
        
        for n=1:20
            if (~isnan(max(ax(n,:))))
                main_ax(n,:)=find(ax(n,:)==max(ax(n,:)));
            else
                main_ax(n,:)=NaN;
            end
        end
        
        for k=1:3
            dum(1,k)=numel(find(main_ax==k));
        end
        
        maxim=numel(find(~isnan(main_ax)));
        if maxim~=max(dum)
            prob(sub,pp)=(max(dum)./maxim);
        else
            prob(sub,pp)=1;
        end
        all(:,pp)=main_ax;
        clear ax main_ax maxim dum
    end
    
    fig=figure(1)
    subplot(2,5,sub)
    degs=0:30:330;
    plot(degs,prob(sub,:),'r.','MarkerSize',10)
    ylabel({'probability of consistent ','main axis'})
    xlabel('stimulation phase (degrees)')
    xlim([-5 335])
    xticks([0:30:330])
    ylim([0 1])
    set(gca,'FontSize',10)
    set(gca,'FontName','Arial','XTickLabelRotation',45)
    set(fig,'color','w');
    box('off')
    
    clearvars -except out prob sub
end


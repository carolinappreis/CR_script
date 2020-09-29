function [out]=cvar(out)
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/seg_envelope.mat')

for iii=1:size(out.start_c,1)
    m_ax=1;
    ARC=squeeze(out.change_c{iii,2}{m_ax,1});
    sup1=[];
    amp1=[];
    
    for yi=1:size(ARC,2)
        
      
        dum=(env_seg{1,yi})';
        dum=(dum(:))';
        if ARC(iii,yi)>0
            amp1=[amp1  dum];
        else
            sup1=[sup1  dum];
        end
    end
    
    %     amp=padarray(amp1,[0 2000],0,'both');
    %     sup=padarray(sup1,[0 2000],0,'both');
    
    clear dum
    
    cv_s(iii,:)=nanstd(sup1)/nanmean(sup1);
    
    cv_a(iii,:)=nanstd(amp1)/nanmean(amp1);
    
    cv_ns(iii,:)=nanstd(env_ns{iii,1})/nanmean(env_ns{iii,1}); 
    
    figure(1)
    subplot(2,5,iii)
    plot(sup1,'LineWidth',2)
    hold on
    plot(amp1,'LineWidth',2)
    plot(env_ns{iii,1},'LineWidth',2)
  
    clear amp1 sup1
    
end

[p,stats.cv_ns_s]=ttest(cv_ns,cv_s);
[p,stats.cv_ns_a]=ttest(cv_ns,cv_a);
[p,stats.cv_s_a]=ttest(cv_s,cv_a);

figure(2)
bar([median(cv_ns) median(cv_s) median(cv_a)],'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.3,'EdgeColor','none');
hold on
plot([cv_ns cv_s cv_a]','.','MarkerSize',15)
plot([cv_ns cv_s cv_a]','MarkerSize',15)
xlim([0 4])
xticklabels({'NS','SUP','AMP'})
box('off')
ylabel('coeffiecient of variation')
set(gca,'FontSize',12)


end


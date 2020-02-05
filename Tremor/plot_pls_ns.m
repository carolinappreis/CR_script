clear all
close all
 load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\reference_pls_ns.mat')
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/reference_pls_ns.mat')

for i=1:10;
f1=figure(1)
subplot(2,5,i)
    plot((time_ns(i,:)),'.','MarkerSize',15,'Color',[0.5 0.5 0.5])
    hold on
    plot((time_ns(i,:)),'LineWidth',1.5,'Color',[0.5 0.5 0.5])
    
    plot((time_pls(i,:)),'.','MarkerSize',15,'Color',[0 0.5 0.5])
    plot((time_pls(i,:)),'LineWidth',1.5,'Color',[0 0.5 0.5])
    
    xlim([0 4])
    xticks([1 2 3])
    xticklabels({'rest','bef stim','end stim'})
    xtickangle(45)
    ylabel({'Tremor severity';'(zscore)'})
    set(gca,'FontSize',12)
    box('off')
    
end

    
f1.Units = 'centimeters';
f1.OuterPosition= [10, 10, 30, 15];
set(f1,'color','w');
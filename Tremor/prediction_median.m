
clear all
close all
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/3_ax_medsplit_0420.mat')
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/pls_change_all.mat')


stim={[8 10];[5 8];[8];[5];[8];[5];[2];[9 12];[9];[2 12]}; % stim phase-phaselocked

amp_thresh=a_prior; clear a_prior;
change_stim=m_change; clear m_change;
amp_pls=ss; clear ss;
match=NaN(10,3);


for i=1:10
    for ii=1:3
        if amp_pls(i)<amp_thresh(i,ii) &  (~isempty(find(arc1(i,ii,stim{i,1})>0)) && change_stim(i,ii)>0 ||...
                ~isempty(find(arc1(i,ii,stim{i,1})<0)) && change_stim(i,ii)<0 )...
                | ( amp_pls(i)>amp_thresh(i,ii)  & ( ~isempty(find(arc2(i,ii,stim{i,1})>0)) && change_stim(i,ii)>0 ||....
                ~isempty(find(arc2(i,ii,stim{i,1})<0)) && change_stim(i,ii)<0))
            match(i,ii)=1;
        end
    end
end


accuracy=sum(sum(~isnan(match)))/numel(match)*100;

% ARC median splot
for hh=1:size(arc1,1)
    for axx=1:3
        f1=figure(hh)
        subplot(1,3,axx)
        y=[squeeze(arc1(hh,axx,:))';squeeze(arc2(hh,axx,:))']';
        b = bar(0:30:330,y,'EdgeColor','none');
        yline(0,'LineWidth',1)
        ylim([-max(max(y)) max(max(y))])
        %     yticks([ -1:0.25:1])
        box('off')
        ylabel({'Change in tremor severity'})
        xlabel({'Stimulation phase (degrees)'})
        
        f1.Units = 'centimeters';
         f1.OuterPosition= [10, 10, 20, 8];
        set(gca,'XTickLabelRotation',45)
            set(gca,'FontSize',8)
        set(f1,'color','w');
    end
    filename=['3ax_medsplit_' num2str(hh) '.pdf'];
%     saveas(gcf,filename)
    
end

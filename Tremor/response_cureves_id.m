
clear all
close all

metric=0;

if metric==0;
            load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/am_ax')
            load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/amp_ARC','LS')
            load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/amp_NS','no_s')
            load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/smooth_arc3')
        load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','squash')
          S=am_ax;
    
%     load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\am_ax')
%     load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\amp_NS','no_s')
%     load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\amp_ARC','LS')
%     load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\smooth_arc3')
%     load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','blushred','squash');
%     S=am_ax;
    
    cl=blushred;
    cl1=squash;
    
else
            load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/fm_ax')
            load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/freq_FRC','LS')
            load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/freq_NS','no_s')
            load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/smooth_frc3')
            load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','aegean','stone');
    
%     load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\fm_ax')
%     load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\freq_FRC','LS')
%     load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\freq_NS','no_s')
%     load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\smooth_frc3')
%     load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','aegean','stone');
    S=fm_ax;
    
    cl=aegean;
    cl1=stone;
end

for nr=1
    close all
    for i=1:size(smo_s,1)
        if nr==1
            f1=figure(i)
            bar((0:30:330),smo_s(i,:),'LineWidth',0.5,'FaceColor',cl,'EdgeColor',cl)
            yline(0,'LineWidth',1)
            box('off')
            ylabel('Change in tremor severity')
            xlabel('Stimulation phase (degrees)')
            f1.Units = 'centimeters';
            f1.OuterPosition= [10, 10, 12, 12];
            set(gca,'XTickLabelRotation',45)
            set(gca,'FontSize',14)
            set(f1,'color','w');
%             cd('C:\Users\creis\Desktop\arcs_share\peripheral')
                        cd('/Users/Carolina/OneDrive - Nexus365/arcs_share/peripheral')

            saveas(['arc_',num2str(i)],'.png')
            
        elseif nr==2
            f1=figure(i)
            bar((0:30:330),smo_s(i,:),'LineWidth',0.5,'FaceColor',cl,'EdgeColor',cl)
            hold on
            bar((0:30:330),smo_s1(i,:),'LineStyle','--','LineWidth',1,'FaceColor','none','EdgeColor','k')
            yline(0,'LineWidth',1)
            box('off')
            ylabel('Change in tremor severity')
            xlabel('Stimulation phase (degrees)')
            f1.Units = 'centimeters';
            f1.OuterPosition= [10, 10, 12, 12];
            set(gca,'XTickLabelRotation',45)
            set(gca,'FontSize',14)
            set(f1,'color','w');
            cd('C:\Users\creis\Desktop\arcs_share\peripheral')
            savefig(['arc_trials_mainaxis_',num2str(i)])
            
        elseif nr==3
            f1=figure(i)
            bar((0:30:330),smo_s(i,:),'FaceColor',cl,'EdgeColor',cl)
            hold on
            plot((0:30:330),smo_s2(i,:),'LineWidth',1,'Color','k')
            plot((0:30:330),smo_s3(i,:),'LineWidth',1,'Color','k')
            yline(0,'LineWidth',1)
            box('off')
% %             ylabel('Change in tremor severity')
            ylabel('Change in instataneous frequency')
            xlabel('Stimulation phase (degrees)')
            f1.Units = 'centimeters';
            f1.OuterPosition= [10, 10, 12, 12];
            set(gca,'XTickLabelRotation',45)
            set(gca,'FontSize',14)
            set(f1,'color','w');
%             cd('C:\Users\creis\Desktop\arcs_share\peripheral')
            cd('/Users/Carolina/OneDrive - Nexus365/arcs_share/peripheral')

            savefig(['frc_plus2axes',num2str(i)])
            
        elseif nr==4
            f1=figure(i)
            bar((0:30:330),S(i,:),'FaceColor',cl,'EdgeColor',cl)
            hold on
%             yline(prctile(no_s(i,:),99.7917),'k--','LineWidth',1)
%             yline(prctile(no_s(i,:),0.2083),'k--','LineWidth',1)
            yline(prctile(LS(i,:),99.7917),'k-.','LineWidth',1)
            yline(prctile(LS(i,:),0.2083),'k-.','LineWidth',1)
            ylabel('Change in tremor severity')
            xlabel('Stimulation phase (degrees)')
            f1.Units = 'centimeters';
            f1.OuterPosition= [10, 10, 12, 12];
            set(gca,'XTickLabelRotation',45)
            set(gca,'FontSize',14)
            set(f1,'color','w');
            box('off')
%             cd('C:\Users\creis\Desktop\arcs_share\peripheral')
%             cd('/Users/Carolina/OneDrive - Nexus365/arcs_share/peripheral')
% 
%             savefig(['arc_raw_randstim_',num2str(i)])
        elseif nr==5
            f1=figure(i)
            bar((0:30:330),smo_s(i,:),'LineWidth',0.5,'FaceColor',cl,'EdgeColor',cl)
            yline(0,'LineWidth',1)
            box('off')
            ylabel('Change in tremor severity')
            xlabel('Stimulation phase (degrees)')
            f1.Units = 'centimeters';
            f1.OuterPosition= [10, 10, 12, 12];
            set(gca,'XTickLabelRotation',45)
            set(gca,'FontSize',14)
            set(f1,'color','w');
            cd('/Users/Carolina/OneDrive - Nexus365/arcs_share/peripheral')
            savefig(['frc_',num2str(i)])
        end
    end
end


% load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','blushred','squash');
%  load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','squash')
clear all
close all
load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\arc_mediansplit.mat')
load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','blushred','squash');

cl=blushred;
cl1=squash;

for i=1:size(arc1,1)
    f1=figure(i)
    bar(arc2(i,:),'FaceColor',cl,'EdgeColor',cl)
    hold on
    bar(arc1(i,:),'LineStyle','--','LineWidth',1,'FaceColor','none','EdgeColor','k')
    ylabel('Change in tremor severity')
    xlabel('Stimulation phase (degrees)')
    f1.Units = 'centimeters';
    f1.OuterPosition= [10, 10, 12, 12];
    box('off')
    set(gca,'XTickLabelRotation',45)
    set(gca,'FontSize',14)
    set(f1,'color','w');
    cd('C:\Users\creis\Desktop\arcs_share\peripheral')
    savefig(['arc_Amp_mediasplit',num2str(i)])
end






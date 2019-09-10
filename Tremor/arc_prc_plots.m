
clear all
close all

metric=1;

if metric==0;
    % load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/smooth_arc3')
    % load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','squash')
    
    load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\am_ax')
    load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\amp_NS','no_s')
    load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\amp_ARC','LS')
    load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\smooth_arc3')
    load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','blushred','squash');
    S=am_ax;
    
    cl=blushred;
    cl1=squash;
    
else
    %     load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/fm_ax')
    %     load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/freq_FRC','LS')
    %     load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/freq_NS','no_s')
    %     load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/smooth_frc3')
    %     load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','aegean','stone');
    %
    load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\fm_ax')
    load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\freq_FRC','LS')
    load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\freq_NS','no_s')
    load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\smooth_frc3')
    load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','aegean','stone');
    S=fm_ax;
    
    cl=aegean;
    cl1=stone;
end

for i=1:size(smo_s,1)
    f1=figure(1)
    subplot(1,size(smo_s,1),i)
    bar(smo_s(i,:),'LineWidth',0.5,'FaceColor',cl,'EdgeColor',cl)
    hold on
    bar(smo_s1(i,:),'LineStyle','--','LineWidth',1,'FaceColor','none','EdgeColor','k')
    yline(0,'LineWidth',1)
    box('off')

    f2=figure(2)
    subplot(1,size(smo_s,1),i)
    bar(smo_s(i,:),'FaceColor',cl,'EdgeColor',cl)
    hold on
    plot(smo_s2(i,:),'LineWidth',1,'Color','k')
    plot(smo_s3(i,:),'LineWidth',1,'Color','k')
    yline(0,'LineWidth',1)
    box('off')

    f3=figure(3)
    subplot(1,size(S,1),i)
    bar(S(i,:),'FaceColor',cl,'EdgeColor',cl)
    hold on
    yline(prctile(no_s(i,:),99.7917),'r--','LineWidth',1)
    yline(prctile(no_s(i,:),0.2083),'r--','LineWidth',1)
    
    yline(prctile(LS(i,:),99.7917),'k--','LineWidth',1)
    yline(prctile(LS(i,:),0.2083),'k--','LineWidth',1)
    
    box('off')
    
end


f1.Units = 'centimeters';
f1.OuterPosition= [10, 10, 60, 6];
set(gca,'FontSize',8)
set(f1,'color','w');

f2.Units = 'centimeters';
f2.OuterPosition= [10, 10, 60, 6];
set(gca,'FontSize',8)
set(f2,'color','w');
%
f3.Units = 'centimeters';
f3.OuterPosition= [10, 10, 60, 6];
set(gca,'FontSize',8)
set(f3,'color','w');


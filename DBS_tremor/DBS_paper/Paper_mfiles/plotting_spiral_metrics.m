

% cd('/Users/Carolina/Library/CloudStorage/OneDrive-Nexus365/Phasic_DBS/patient data/DBS_DATA/clean_SH_spirals/SPSS')
% 
% cd('/Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/DBS_paper')

clear; close all

load('/Users/Carolina/Library/CloudStorage/OneDrive-Nexus365/Phasic_DBS/patient data/DBS_DATA/clean_SH_spirals/SPSS/group_acc_means.mat')
load('/Users/Carolina/Library/CloudStorage/OneDrive-Nexus365/Phasic_DBS/patient data/DBS_DATA/clean_SH_spirals/SPSS/group_t_means.mat')



 load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','jam','yam','teal','blushred','aegean','stone','squash','sapphire','azure','purple','berryred', 'rosewood');
% color_b1=[aegean;blushred;sapphire];
 color_b1=[jam;teal;yam];

% CT=cbrewer('div','RdGy',9);
% 
% color_b1=[CT(1,:); CT(2,:); CT(3,:)];


for i=1:2
    if i==1
        dat=stim;
    else
        dat=stim_t;
    end
    err=[dat(1,2) dat(3,2) dat(2,2)];
    x=[dat(1,1) dat(3,1) dat(2,1)];
    plot([1:3],x,'--','Color',color_b1(i,:),'LineWidth',1)
    hold on
    errorbar([1:3],x,err,'o','Color',color_b1(i,:),'LineWidth',1.5,'MarkerSize',5,'MarkerEdgeColor',color_b1(i,:),'MarkerFaceColor',color_b1(i,:))
    xlim([0.75 3.25])
    ylim([0 2.5])
    yticks(0:0.5:2)
    xticks(1:3)
    xticklabels({'NS','PLS','HFS'})
    xlabel ('Stimualtion conditions')
    ylabel('Mean tremor envelope')
    box('off')
    set(gca,'FontSize',14)
    clear dat
end


for uu=1:3
cd('/Users/Carolina/Library/CloudStorage/OneDrive-Nexus365/Phasic_DBS/patient data/DBS_DATA/clean_SH_spirals/SPSS')
fname=strcat('sq_',(num2str(uu)),'.xlsm');
sq=readtable(fname);
figure(2)
subplot(1,3,uu)
dat= [sq.Mean(1:3:10) sq.Mean((1:3:10)+1) sq.Mean((1:3:10)+2)];
d_err=[sq.Std_Error(1:3:10) sq.Std_Error((1:3:10)+1) sq.Std_Error((1:3:10)+2)];
for i=1:3
    x=dat(:,i); err=d_err(:,i);
    plot([1:4],x,'-','Color',color_b1(i,:),'LineWidth',1.5)
    hold on
    errorbar([1:4],x,err,'o','Color',color_b1(i,:),'LineWidth',1,'MarkerSize',3,'MarkerEdgeColor',color_b1(i,:),'MarkerFaceColor',color_b1(i,:))
     xlim([0.75 4.25])
      ylim([0 2.9])
    yticks(0:0.5:2.5)
     xticks(1:4)
    xticklabels({'Q1','Q2','Q3','Q4'})
%     xlabel ('Spiral polar quadrant')
%     ylabel('Mean tremor envelope')
    box('off')
    title(strcat('patient',num2str(uu+1)))
    set(gca,'FontSize',12)
    clear x err
end
clear fname sq dat d_err
end
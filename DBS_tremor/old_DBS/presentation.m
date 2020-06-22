
clear all
close all

        load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/am_ax')
        load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_amp_ARC','LS','tt1','ttall')        
        load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/NS_PS_result','idv_NS')
        load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/smooth_arc3')
        load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','squash');

    
    S=ttall;
    cl=blushred;
    cl1=squash;
    
% close all
% i=1
%     f1=figure;
%     y= smo_s(i,:);
%     
%     f1=figure(1)
%     bar(y,'LineWidth',0.5,'FaceColor',cl,'EdgeColor',cl)
%     hold on
%     xticks([]);
%     yticks([]);
%     hold on
%     rsg=sin_fit1(y);
%     rs_gauss(i,:)=rsg.adjrsquare;
%     xlabel('');ylabel('');legend('off')
%     ylabel('Change in tremor severity')
%     xlabel('Stimulation phase (degrees)')
%     f1.Units = 'centimeters';
%     f1.OuterPosition= [10, 10, 12, 12];
%     box('off')
%     set(gca,'XTickLabelRotation',45)
%     set(gca,'FontSize',14)
%     set(f1,'color','w');   
%     

close all
i=1
cl=[0.5 0.5 0.5];
f1=figure;
y= tt1{1, 1}{1, 1};
f1=figure(1)
bar(0:30:330,100.*nanmedian(y),'LineWidth',0.5,'FaceColor',cl,'EdgeColor',cl)
hold on
plot(0:30:330,100.*y,'.')
ylabel('Change in tremor severity')
xlabel('Stimulation phase (degrees)')
f1.Units = 'centimeters';
f1.OuterPosition= [10, 10, 12, 12];
set(gca,'XTickLabelRotation',45)
box('off')
set(gca,'FontSize',14)
set(f1,'color','w');


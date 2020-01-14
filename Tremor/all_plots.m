

clear all
close all
load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\cleaned_rc12.mat')
load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\newnonstim2.mat')
load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\labels_shift.mat')
load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\arc_mediansplit_0120.mat')
load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','blushred','squash');
% 
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cleaned_rc12.mat')
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/newnonstim2.mat')
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/labels_shift.mat')
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/arc_mediansplit_0120.mat')
% load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','squash');
cl=blushred;

%attention to number 13 and 14 of tt1 (they are the same pateint, 13 5
%trials with 5 pulses stim and 14, 5 trials with 1 pulse stim


%smoothed ARC'S
% iiii=[1 2 3 4 5 8 10 11 12 13 16 17 19 20];

cohort=[2 3 4 5 8 10 11 13 16 17];
dum=intersect(iiii,cohort);

pt=[];
for i=1:length(dum)
    pt=[pt find(iiii==dum(i))];
end

  main=[1 1 3 1 3 3 3 3 1 1];
%   main=[1 1 1 1 1 1 1 1 1 1];
for pp=1:size(pt,2)
    for kk=1:size(tt1,2)
        
        curve=nanmedian(tt1{pt(pp),kk},1);
        SO=repmat(curve,1,3);
        
        for i=size(tt1{pt(pp),kk},2)+1:size(tt1{pt(pp),kk},2)*2
            smooth_c(1,i-12)=sum(SO(1,(i-1:i+1)))./length(SO(1,(i-1:i+1)));
        end
        
        
        smoo_all(pp,kk,:)=smooth_c;  clear smooth_c
        nostim_all(pp,kk,:)=nostim(pp,kk,:);
        nostim_var(pp,kk,:)=var(nostim(pp,kk,:));
        test(pp,kk,:)=prctile(nostim(pp,kk,:),99.7917);
        raw(pp,kk,:)=nanmedian(tt1{pt(pp),kk});
        ray_p(pp,kk,1)=circ_rtest(squeeze(nanmedian(tt1{pt(pp),kk})));
        ray_test(pp,kk,1)=(circ_rtest(squeeze(nanmedian(tt1{pt(pp),kk}))))<0.05;
    
    end
   
%     bar(nostim_var)
%     figure(2)
%     bar(test)
% %     dum=(squeeze(nostim_all(pp,:,:)))';
%     boxplot(dum,'PlotStyle','compact')
    
    
     
    smoo_main(pp,:)=squeeze(smoo_all(pp,1,:));
    raw_main(pp,:)=squeeze(nanmedian(tt1{pt(pp),1}));
    nostim_1(pp,:)=squeeze(nostim(pp,main(pp),:));
    label_shift_1(pp,:)=squeeze(LS(pp,main(pp),:));
end



for i =1:size(smoo_main,1)
    
%     %smooth arc main axes
%     f1= figure(1) 
%     subplot(2,5,i)
%     bar(0:30:330,smoo_main(i,:),'FaceColor',cl,'EdgeColor',cl)
%     ylabel('Change in tremor severity')
%     xlabel('Stimulation phase (degrees)')
%     set(gca,'XTickLabelRotation',45)
%     set(gca,'FontSize',12)
%     box('off')
%     
%     
%     %raw ARC with NS thereshold for main axes
%     f2= figure(2)
%     subplot(2,5,i)
%     bar(0:30:330,raw_main(i,:),'FaceColor',cl,'EdgeColor',cl)
%     hold on
%     yline(prctile(nostim_1(i,:),99.7917),'k--','LineWidth',1)
%     yline(prctile(nostim_1(i,:),0.2083),'k--','LineWidth',1)
%     ylabel('Change in tremor severity')
%     xlabel('Stimulation phase (degrees)')
%     set(gca,'XTickLabelRotation',45)
%     set(gca,'FontSize',12)
%     box('off')
%     
%     
%     f3=figure(3);
%     subplot(2,5,i)
%     bar(0:30:330,squeeze(smoo_all(i,1,:)),'FaceColor',cl,'FaceAlpha',[0.3],'EdgeColor','none')
%     hold on 
%     plot(0:30:330,squeeze(smoo_all(i,2,:)),'k','LineWidth',1.5)
%     plot(0:30:330,squeeze(smoo_all(i,3,:)),'k','LineWidth',1.5)
%     ylabel('Change in tremor severity')
%     xlabel('Stimulation phase (degrees)')
%     set(gca,'XTickLabelRotation',45)
%     box('off')
    
    %raw ARC with Phase shift thereshold for main axes
    f4= figure(4)
    subplot(2,5,i)
    bar(0:30:330,raw_main(i,:),'FaceColor',cl,'EdgeColor',cl)
    hold on
    yline(prctile(label_shift_1(i,:),99.7917),'k--','LineWidth',1)
    yline(prctile(label_shift_1(i,:),0.2083),'k--','LineWidth',1)
    ylabel('Change in tremor severity')
    xlabel('Stimulation phase (degrees)')
    set(gca,'XTickLabelRotation',45)
    set(gca,'FontSize',12)
    box('off')
    
end

f1.Units = 'centimeters';f2.Units = 'centimeters';f3.Units = 'centimeters';f4.Units = 'centimeters';
f1.OuterPosition= [10, 10, 55, 15]; f2.OuterPosition= [10, 10, 50, 15];f3.OuterPosition= [10, 10, 50, 15];f4.OuterPosition= [10, 10, 50, 15];
set(f1,'color','w'); set(f2,'color','w');set(f3,'color','w');set(f4,'color','w');



%raw ARC with NS thereshold for 3 axis

for i=1:size(pt,2)
    f1=figure()
    for axis=1:3
        subplot(1,3,axis)
        bar(0:30:330,squeeze(raw(i,axis,:)),'FaceColor',cl,'EdgeColor',cl)
        hold on
        yline(prctile(nostim_1(i,:),99.7917),'k--','LineWidth',1)
        yline(prctile(nostim_1(i,:),0.2083),'k--','LineWidth',1)
        
        ylabel('Change in tremor severity')
        xlabel('Stimulation phase (degrees)')
        f1.Units = 'centimeters';
        f1.OuterPosition= [10, 10, 12, 12];
        set(gca,'XTickLabelRotation',45)
        box('off')
           
    end
    f1.Units = 'centimeters';
    f1.OuterPosition= [10, 10, 30, 10];
    set(f1,'color','w');
      
end



% ARC median splot
for hh=1:size(arc1,1)
    f1=figure(1)
    subplot(2,5,hh)
    y=[arc1(hh,:);arc2(hh,:)]';
    b = bar(0:30:330,y);
    yline(0,'LineWidth',1)
%     ylim([-0.75 0.75])
%     yticks([ -1:0.25:1])
    box('off')
    ylabel({'Change in tremor severity'})
    xlabel({'Stimulation phase (degrees)'})
    
    f1.Units = 'centimeters';
    f1.OuterPosition= [10, 10, 50, 15];
    set(gca,'XTickLabelRotation',45)
    set(gca,'FontSize',12)
    set(f1,'color','w');   
end




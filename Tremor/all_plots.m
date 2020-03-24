

clear all
close all
% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\cleaned_rc12_noaddon.mat')
% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\newnonstim10.mat')
% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\labels_shift.mat','LS')
% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\ns_median_split.mat','ns_msplit')
% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\arc_mediansplit_0220.mat','arc1','arc2')
% 
% load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','blushred','squash');

load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cleaned_rc12_noaddon.mat')
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/newnonstim10.mat')
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/labels_shift.mat')
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/ns_median_split.mat')
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/arc_mediansplit_0220.mat')
load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','squash');
cl=blushred;

%attention to number 13 and 14 of tt1 (they are the same pateint, 13 5
%trials with 5 pulses stim and 14, 5 trials with 1 pulse stim


% smoothed ARC'S
% iiii=[1 2 3 4 5 8 10 11 12 13 16 17 19 20];

% cohort=[2 3 4 5 8 10 11 13 16 17];
% dum=intersect(iii,cohort);
%
% pt=[];
% for i=1:length(dum)
%     pt=[pt find(iii==dum(i))];
% end
%
% for i=1:size(pt,2)
%     for ii=1:3
%     ns(i,ii,:)=nostim(pt(i),ii,:);
%     nsout(i,ii,:)=nostimout(pt(i),ii,:);
%     end
% end
%
% clear nostim
% nostim=ns;
% clear nostimout
% nostimout=nsout;


main=[1 1 3 1 3 3 3 3 1 1];
ns_mat={[1 2 3];[1 2 3];[3 2 1];[1 2 3];[3 2 1];[3 2 1];[3 2 1];[3 2 1];[1 2 3];[1 2 3]};

for ik=1:size(tt1,1)
    for kk=1:size(tt1,2)
        
        curve=nanmedian(tt1{ik,kk},1);
        SO=repmat(curve,1,3);
        
        for i=size(tt1{ik,kk},2)+1:size(tt1{ik,kk},2)*2
            smooth_c(1,i-12)=sum(SO(1,(i-1:i+1)))./length(SO(1,(i-1:i+1)));
        end
        
        
        smoo_all(ik,kk,:)=smooth_c;  clear smooth_c
        nostim_var(ik,kk,:)=var(nostim(ik,kk,:));
        raw(ik,kk,:)=nanmedian(tt1{ik,kk});
        
        %         p (ik,kk,1) = kruskalwallis(tt1{ik,kk});
        ns_mat1(ik,kk,:)=squeeze(nostim(ik,ns_mat{ik,1}(kk),:));
    end
    
    %     bar(nostim_var)
    %     figure(2)
    %     bar(test)
    % %     dum=(squeeze(nostim_all(ik,:,:)))';
    %     boxplot(dum,'PlotStyle','compact')
    
    smoo_main(ik,:)=squeeze(smoo_all(ik,1,:));
    raw_main(ik,:)=squeeze(nanmedian(tt1{ik,1}));
    nostim_1(ik,:)=squeeze(nostim(ik,main(ik),:));
    label_shift_1(ik,:)=squeeze(LS(ik,main(ik),:));
    
end

mean(nanmean(vertcat(tt1{1:10,1})))*100
std(nanmean(vertcat(tt1{1:10,1})))*100

for i =1:size(smoo_main,1)
    
    % %     %smooth arc main axes
    f1= figure(1)
    subplot(2,5,i)
    bar(0:30:330,smoo_main(i,:),'FaceColor',cl,'EdgeColor',cl)
    ylabel('Change in tremor severity')
    xlabel('Stimulation phase (degrees)')
    set(gca,'XTickLabelRotation',45)
    box('off')
    %
    %
    %     %raw ARC with NS thereshold for main axes
    f2= figure(2)
    subplot(2,5,i)
    bar(0:30:330,raw_main(i,:),'FaceColor',cl,'EdgeColor',cl)
    hold on
    yline(prctile(nostim_1(i,:),99.7917),'k--','LineWidth',1)
    yline(prctile(nostim_1(i,:),0.2083),'k--','LineWidth',1)
    ylabel('Change in tremor severity')
    xlabel('Stimulation phase (degrees)')
    ylim([-0.8 0.8])
    set(gca,'XTickLabelRotation',45)
    box('off')
    %
    
    f3=figure(3);
    subplot(2,5,i)
    bar(0:30:330,squeeze(smoo_all(i,1,:)),'FaceColor',cl,'FaceAlpha',[0.3],'EdgeColor','none')
    hold on
    plot(0:30:330,squeeze(smoo_all(i,2,:)),'k','LineWidth',1.5)
    plot(0:30:330,squeeze(smoo_all(i,3,:)),'k','LineWidth',1.5)
    ylim([-0.8 0.8])
    ylabel('Change in tremor severity')
    xlabel('Stimulation phase (degrees)')
    set(gca,'XTickLabelRotation',45)
    box('off')
    
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
    box('off')
    
    up(i,1)=prctile(nostim_1(i,:),99.7917);
    down(i,1)=prctile(nostim_1(i,:),0.2083);
end

f1.Units = 'centimeters';f2.Units = 'centimeters';f3.Units = 'centimeters';f4.Units = 'centimeters';
f1.OuterPosition= [10, 10, 55, 15]; f2.OuterPosition= [10, 10, 50, 15];f3.OuterPosition= [10, 10, 50, 15];f4.OuterPosition= [10, 10, 50, 15];
set(f1,'color','w'); set(f2,'color','w');set(f3,'color','w');set(f4,'color','w');





% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\cleaned_rc12_noaddon.mat')
% load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','blushred','squash');
% cl=blushred;
% for i=1:10
%     figure(1)
%     subplot(2,5,i)
%     bar(0:30:330,raw_main(i,:),'FaceColor',cl,'EdgeColor',cl)
%     hold on
%     bar(0:30:330,nanmedian(tt1{i,1}),'FaceColor','none','EdgeColor','k')
%     ylabel('Change in tremor severity')
%     xlabel('Stimulation phase (degrees)')
%     set(gca,'XTickLabelRotation',45)
%     set(gca,'FontSize',12)
%     box('off')
% end
%





%raw ARC with NS thereshold for 3 axis

for i=1:size(tt1,1)
    f1=figure()
    for axis=1:size(tt1,2)
        subplot(1,3,axis)
        bar(0:30:330,squeeze(raw(i,axis,:)),'FaceColor',cl,'EdgeColor',cl)
        hold on
        yline(prctile(ns_mat1(i,axis,:),99.7917),'k--','LineWidth',1)
        yline(prctile(ns_mat1(i,axis,:),0.2083),'k--','LineWidth',1)
        
        ylabel('Change in tremor severity')
        xlabel('Stimulation phase (degrees)')
        f1.Units = 'centimeters';
        f1.OuterPosition= [10, 10, 12, 12];
        set(gca,'XTickLabelRotation',45)
        box('off')
        
    end
    f1.Units = 'centimeters';
    f1.OuterPosition= [10, 10, 30, 8];
    set(f1,'color','w');
    
end


% ARC median splot
for hh=1:size(arc1,1)
    f1=figure(1)
    subplot(2,5,hh)
    y=[arc1(hh,:);arc2(hh,:)]';
    b = bar(0:30:330,y);
    yline(0,'LineWidth',1)
        ylim([-0.75 0.75])
    %     yticks([ -1:0.25:1])
    box('off')
    ylabel({'Change in tremor severity'})
    xlabel({'Stimulation phase (degrees)'})
    
    f1.Units = 'centimeters';
    f1.OuterPosition= [10, 10, 50, 15];
    set(gca,'XTickLabelRotation',45)
%     set(gca,'FontSize',12)
    set(f1,'color','w');
end


%%% ARC median splot plus thresh

for hh=1:size(arc1,1)
    f1=figure(1)
    subplot(2,5,hh)
    y=arc1(hh,:);
    bar(0:30:330,y);
    hold on
    yline(prctile(ns_msplit(hh,1,:),99.7917),'k--','LineWidth',1)
    yline(prctile(ns_msplit(hh,1,:),0.2083),'k--','LineWidth',1)
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
    
    
    if max(y)>= prctile(ns_msplit(hh,1,:),99.7917) | min(y)<=prctile(ns_msplit(hh,1,:),0.2083)
        sig_a1(hh,1,:)=1;
    else
        sig_a1(hh,1,:)=0;
    end
    
    
    f2=figure(2)
    subplot(2,5,hh)
    y=arc2(hh,:);
    bar(0:30:330,y);
    hold on
    yline(prctile(ns_msplit(hh,2,:),99.7917),'k--','LineWidth',1)
    yline(prctile(ns_msplit(hh,2,:),0.2083),'k--','LineWidth',1)
    yline(0,'LineWidth',1)
    %     ylim([-0.75 0.75])
    %     yticks([ -1:0.25:1])
    box('off')
    ylabel({'Change in tremor severity'})
    xlabel({'Stimulation phase (degrees)'})
    
    f2.Units = 'centimeters';
    f2.OuterPosition= [10, 10, 50, 15];
    set(gca,'XTickLabelRotation',45)
    set(gca,'FontSize',12)
    set(f2,'color','w');
    
        if max(y)>= prctile(ns_msplit(hh,2,:),99.7917) | min(y)<=prctile(ns_msplit(hh,2,:),0.2083)
        sig_a2(hh,1,:)=1;
    else
        sig_a2(hh,1,:)=0;
    end
    
end


% ARC all trials plots
for hh=1:size(tt1,1)
    f1=figure(1)
    subplot(2,5,hh)
    bar(0:30:330,nanmedian(tt1{hh,1}),'FaceColor',cl,'EdgeColor',cl);
    hold on
    plot(0:30:330,(tt1{hh,1}),'k.');
    box('off')
    ylabel({'Change in tremor severity'})
    xlabel({'Stimulation phase (degrees)'})
    
    f1.Units = 'centimeters';
    f1.OuterPosition= [10, 10, 50, 15];
    set(gca,'XTickLabelRotation',45)
    set(f1,'color','w');
end


% cr=squeeze(tt1{2,1});
% cr1(1:10,1:12)=cr(1:10,1:12);
% bins=[min(min(cr)) nanmedian(nanmedian(cr)) max(max(cr))];
% [a1,ef1]=find(cr1<bins(2));
% [a2,ef2]=find(cr1>=bins(2));
% subplot(1,2,1)
% rose(ang(ef1));
% subplot(1,2,2)
% rose(ang(ef2));
% [circ_rtest(ang(ef1)); circ_rtest(ang(ef2))]
% figure;
% bar(0:30:330,nanmedian(cr))
% hold on
% plot(0:30:330,cr,'.')
% yline(bins(2))


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
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/labels_shift.mat')
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/ns_median_split.mat')
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/arc_mediansplit_0220.mat')
% load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','squash');
% cl=blushred;

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

d=double(raw_main);
clear nostim; nostim=nostim_1;
n=[];
for i=1:size(tt1,1)
    %     f1 = figure(i)       
        data=d(i,:);
%         n=[n ; data(find(data > prctile(nostim(i, :), 99.7917) | data< prctile(nostim(i, :), 0.2083)))];
        n_phases(i)=numel(find(data > prctile(nostim(i, :), 99.7917) | data< prctile(nostim(i, :), 0.2083)));
        if ~isempty(find(data > prctile(nostim(i, :), 99.7917) | data< prctile(nostim(i, :), 0.2083)))
            phase{i}=find(data > prctile(nostim(i,:), 99.7917) | data< prctile(nostim(i, :), 0.2083));
        else
            phase{i}=NaN;
        end
        
%         subplot(1, 3, axis)
%         bar(0:30:330, squeeze(tt1(i, axis,:)),'FaceColor',cl,'EdgeColor',cl)
%         hold on
%         yline(prctile(nostim(i, axis, :), 99.7917),'k--','LineWidth',1)
%         yline(prctile(nostim(i, axis, :), 0.2083),'k--','LineWidth',1)
%         ylim([-1 1])
%         ylabel('Change in tremor severity')
%         xlabel('Stimulation phase (degrees)')
%         set(gca,'XTickLabelRotation',45)
%         box('off')
        
    
    
    f2= figure(20)
    load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','stone','squash');
    color_b1=[blushred; aegean; stone];
    p=2;
    %ARC sem/std patch
    subplot(2,5,i)
    dr=squeeze(tt1{i,1});
    y2=nanmedian(dr,1);
    y1=prctile(dr,75);
    y3=prctile(dr,25);
    time=0:30:330;
    patch([time fliplr(time)], [y1 fliplr(y2)],[color_b1(p,:)],'FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
    hold on
    patch([time fliplr(time)], [y2 fliplr(y3)],[color_b1(p,:)],'FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
    % plot(y2,'.', 'MarkerSize',20,'Color',color_b1(p,:))
    stem(time,y2,'.', 'LineWidth',1,'MarkerSize',10,'Color',color_b1(p,:))
    yline(0)
    if ~isnan(phase{1,i})
    plot(time(phase{1,i}),y2(phase{1,i}),'.','Color',color_b1(1,:),'MarkerSize',10)
    end
    ylim([-(max(y1)) max(y1)])
    xlim([-5 335])
    xticks([0:30:330])
    box('off')
    ylabel({'Change in tremor severity'})
    xlabel({'Stimulation phase (degrees)'})
    set(gca,'FontSize',12)
    set(gca,'FontName','Arial','XTickLabelRotation',45)
    
    
    %     f1.Units = 'centimeters';
    %     f1.OuterPosition= [10, 10, 30, 8];
    %     set(f1,'color','w');
    f2.Units = 'centimeters';
    f2.OuterPosition= [0, 0, 40, 15];
    set(f2,'color','w');
end







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
    b = bar(0:30:330,y,'EdgeColor','none');
    yline(0,'LineWidth',1)
     ylim([-max(max(y)) max(max(y))])
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



load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','stone','squash');

f1=figure(1);
color_b1=[blushred; aegean; stone];
p=2
for hh=1:size(tt1,1)
    subplot(2,5,hh)
dr=tt1{hh,1};
y2=nanmedian(dr,1);
% sem = nanstd(dr)./ sqrt(size(y2,2));
% stdev = nanstd(dr,1);
stdev=prctile(dr,75);
bar(0:30:330,y2,'FaceColor',color_b1(p,:),'EdgeColor',color_b1(p,:),'FaceAlpha',[0.50],'EdgeAlpha',[0.50]);
hold on
bin=[0:30:330];
for i = 1:length(bin)
    errorbar(bin(i), y2(:,i), stdev(:,i), 'Color',color_b1(p,:), 'linestyle', 'none');
end
ylim([-(max(abs(stdev)+y2)) max(abs(stdev)+y2)])
% ylim([-1.5 1.5])
box('off')
    ylabel({'Change in tremor severity'})
    xlabel({'Stimulation phase (degrees)'})
    f1.Units = 'centimeters';
    f1.OuterPosition= [0, 0, 40, 15];
    set(gca,'FontName','Arial','XTickLabelRotation',45)
    set(f1,'color','w');
    
end



p=2;
f1=figure(1)
%ARC sem/std patch
for hh=1:size(tt1,1)
    subplot(2,5,hh)
    
dr=tt1{hh,1};
y2=nanmedian(dr,1);
y1=prctile(dr,75);
y3=prctile(dr,25);
time=0:30:330;
patch([time fliplr(time)], [y1 fliplr(y2)],[color_b1(p,:)],'FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
hold on
patch([time fliplr(time)], [y2 fliplr(y3)],[color_b1(p,:)],'FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
% plot(y2,'.', 'MarkerSize',20,'Color',color_b1(p,:))
stem(time,y2,'.', 'LineWidth',1,'MarkerSize',10,'Color',color_b1(p,:))
hold on
yline(0)


ylim([-(max(y1)) max(y1)])
xlim([-5 335])
xticks([0:30:330])
% for i = 1:12
%     errorbar(i, y2(:,i), y1(:,i), 'Color',color_b1(p,:), 'linestyle', 'none');
% end
% ylim([-(round(max(stdev+y2),2)) round(max(stdev+y2),2)])
% box('off')

set(gca,'FontSize',12)
box('off')

box('off')
    ylabel({'Change in tremor severity'})
    xlabel({'Stimulation phase (degrees)'})
    f1.Units = 'centimeters';
    f1.OuterPosition= [0, 0, 40, 15];
    set(gca,'FontName','Arial','XTickLabelRotation',45)
    set(f1,'color','w');
end




y2 = nanmedian(dr);
stdev = nanstd(dr)./ sqrt(size(y2,2));
% stem(y2,'.', 'LineWidth',2,'MarkerSize',20,'Color',color_b1{1,1})
plot(y2,'.', 'MarkerSize',20,'Color',color_b1{1,1})
hold on
yline(0)
for i = 1:12
    errorbar(i, y2(:,i), stdev(:,i), 'Color',[0.5 0.5 0.5], 'linestyle', 'none');
end
ylim([-(round(max(stdev+y2),2)) round(max(stdev+y2),2)])
box('off')

ylim([-1.5 1.5])
xlim([1 12])





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
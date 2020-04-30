clear all
close all

method = 'ward-pca2'; %% Cluster method - change it here

c = 2; %% Cluster number - change it here


load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/clusters_BA.mat','nostim_c1')
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/clusters_BA.mat','nostim_c2')
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/output_auto.mat','TT1_C1')
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/output_auto.mat','TT1_C2')


cl = 'r';

%%

if c == 1
    % Cluster 1
    tt1 = TT1_C1;
    nostim = nostim_c1;
else % c == 2
    % Cluster 2
    tt1 = TT1_C2;
    nostim = nostim_c2;
end

cohort = [ 2 3 4 5 8 10 11 13 16 17];

idx = [];
for i = 1:length(cohort)
   nostim_1(i,:) = squeeze(nostim(i, 1, :));
   if isempty(tt1{i,1}) || isnan(nostim_1(i, 1))
       idx = [idx, i];
   end
end
clear i nostim_1

tt1(idx,:) = [];
cohort(idx) = [];
nostim(idx,:,:) = [];

%%%%%%%% INACABADO - MUDAR DE ACORDO COM OS RESULTADOS!
if c == 1
    if strcmp(method, 'ward-pca2')     
        main = [1 1 2 1 1 2 1 3 3 1 1 1];
        ns_mat = {[1 2 3]; [1 2 3]; [2 1 3]; [1 2 3]; [1 2 3]; [2 1 3]; [1 2 3]; [3 2 1]; [3 2 1]; [1 2 3]; [1 2 3]; [1 2 3]};  
%     elseif strcmp(method, 'ward-power')     
%         main = [1 1 1 3 2 1 3 3 3 1];
%         ns_mat = {[1 2 3]; [1 2 3]; [1 2 3]; [3 2 1]; [2 1 3]; [1 2 3]; [3 2 1]; [3 2 1]; [3 2 1]; [1 2 3]};  
%     elseif strcmp(method, 'ward')     
%         main = [1 1 1 1 3 3 1 3 1 1];
%         ns_mat = {[1 2 3]; [1 2 3]; [1 2 3]; [1 2 3]; [3 2 1]; [3 2 1]; [1 2 3]; [3 2 1]; [1 2 3]; [1 2 3]};     
%     elseif strcmp(method, 'average')        
%         main = [1 1 1 3];
%         ns_mat = {[1 2 3]; [1 2 3]; [1 2 3]; [3 2 1]};   
%     elseif strcmp(method, 'complete') 
%         main = [1 1 3 1];
%         ns_mat = {[1 2 3]; [1 2 3]; [3 2 1]; [1 2 3]};
%     elseif strcmp(method, 'single')    
%         main = [1 1 3];
%         ns_mat = {[1 2 3]; [1 2 3]; [3 2 1]};
%     else %% strcmp(method, 'method') 
%         main = [1 1 1 3];
%         ns_mat = {[1 2 3]; [1 2 3]; [1 2 3]; [3 2 1]};
    end
end

if c == 2 %% MUST ADD CASE TO REMOVE FAULTY BASELINE
    if strcmp(method, 'ward-pca2')     
        main = [1 1 1 1 3 2 3 3 1 3 3 1];
        ns_mat = {[1 2 3]; [1 2 3]; [1 2 3]; [1 2 3]; [3 2 1]; [2 1 3]; [3 2 1]; [3 2 1]; [1 2 3]; [3 2 1]; [3 2 1]; [1 2 3]};  
%     elseif strcmp(method, 'ward-power')     
%         main = [1 1 2 1 3 1 3 1 1];
%         ns_mat = {[1 2 3]; [1 2 3]; [2 1 3]; [1 2 3]; [3 2 1]; [1 2 3]; [3 2 1]; [1 2 3]; [1 2 3]};  
%     elseif strcmp(method, 'ward')     
%         main = [];
%         ns_mat = {[1 2 3]; [2 1 3]; [1 2 3]; [3 2 1]; [2 1 3]; [1 2 3]; [1 2 3]; [3 2 1]; [1 2 3]; [1 2 3]; [1 2 3]};     
%     elseif strcmp(method, 'average')        
%         main = [1 1 2 3 1 1 3 3 1 1 1];
%         ns_mat = {[1 2 3]; [1 2 3]; [2 1 3]; [3 2 1]; [1 2 3]; [1 2 3]; [3 2 1]; [3 2 1]; [1 2 3]; [1 2 3]; [1 2 3]};   
%     elseif strcmp(method, 'complete') 
%         main = [1 1 2 1 3 1 1 1 3 3 1 1];
%         ns_mat = {[1 2 3]; [1 2 3]; [2 1 3]; [1 2 3]; [3 2 1]; [1 2 3]; [1 2 3]; [1 2 3]; [3 2 1]; [3 2 1]; [1 2 3]; [1 2 3]};
%     elseif strcmp(method, 'single')    
%         main = [1 1 2 3 1 3 3 1 1];
%         ns_mat = {[1 2 3]; [1 2 3]; [2 1 3]; [3 2 1]; [1 2 3]; [3 2 1]; [3 2 1]; [1 2 3]; [1 2 3]};
%     else %% strcmp(method, 'method') 
%         main = [1 1 2 1 3 1 3 3 3 1 1 1];
%         ns_mat = {[1 2 3]; [1 2 3]; [2 1 3]; [1 2 3]; [3 2 1]; [1 2 3]; [3 2 1]; [3 2 1]; [3 2 1]; [1 2 3]; [1 2 3]; [1 2 3]};
    end
end
%% Plots

for ik = 1:size(tt1,1)
    for kk = 1:size(tt1,2)
        
        curve = nanmedian(tt1{ik,kk},1);
        SO = repmat(curve,1,3);
        
        for i = size(tt1{ik,kk},2)+1:size(tt1{ik,kk},2)*2
            smooth_c(1,i-12) = sum(SO(1,(i-1:i+1)))./length(SO(1,(i-1:i+1)));
        end
        
        
        smoo_all(ik,kk,:) = smooth_c;  clear smooth_c
        nostim_var(ik,kk,:) = var(nostim(ik,kk,:));
        raw(ik,kk,:) = nanmedian(tt1{ik,kk});
        
        %         p (ik,kk,1) = kruskalwallis(tt1{ik,kk});
        ns_mat1(ik,kk,:) = squeeze(nostim(ik, ns_mat{ik,1}(kk), :));
    end
    
    %     bar(nostim_var)
    %     figure(2)
    %     bar(test)
    % %     dum=(squeeze(nostim_all(ik,:,:)))';
    %     boxplot(dum,'PlotStyle','compact')
    
    smoo_main(ik,:) = squeeze(smoo_all(ik, 1, :));
    raw_main(ik,:) = squeeze(nanmedian(tt1{ik, 1}));
    nostim_1(ik,:) = squeeze(nostim(ik,main(ik), :));
%     label_shift_1(ik,:)=squeeze(LS(ik,main(ik),:));
    
end
% cr = reshape(nostim_1, size(nostim_all, 1)*size(nostim_all, 3), 1); 
% cr = abs(cr');
% cr2 = abs(raw_main_s');
% [p,h] = ranksum(cr2, cr) % dist is not normal by hist and samples unpaired (dif size) so ranksum Mann-Withney 
% 
% cr1 = cr(randi(length(cr), 1, length(cr2)));
% [p,h] = signrank(cr2, cr1)

for i = 1:size(smoo_main, 1)
    
    %smooth arc main axes
    f1 = figure(1) 
    subplot(3, 5, i)
    bar(0:30:330, smoo_main(i, :),'FaceColor',cl,'EdgeColor',cl)
    ylabel('Change in tremor severity')
    xlabel('Stimulation phase (degrees)')
    set(gca, 'XTickLabelRotation', 45)
    set(gca, 'FontSize', 12)
    box('off')
%     
%     
    %raw ARC with NS thereshold for main axes
    f2 = figure(2)
    subplot(3, 5, i)
    bar(0:30:330, raw_main(i, :), 'FaceColor', cl, 'EdgeColor', cl)
    hold on
    yline(prctile(nostim_1(i, :), 99.7917), 'k--', 'LineWidth', 1)
    yline(prctile(nostim_1(i, :), 0.2083), 'k--', 'LineWidth', 1)
    ylabel('Change in tremor severity')
    xlabel('Stimulation phase (degrees)')
    set(gca, 'XTickLabelRotation', 45)
    set(gca, 'FontSize', 12)
    box('off')
%     
%     
    f3 = figure(3);
    subplot(3, 5, i)
    bar(0:30:330, squeeze(smoo_all(i, 1, :)), 'FaceColor', cl, 'FaceAlpha', [0.3], 'EdgeColor', 'none')
    hold on 
    plot(0:30:330, squeeze(smoo_all(i, 2, :)), 'k', 'LineWidth', 1.5)
    plot(0:30:330, squeeze(smoo_all(i, 3, :)), 'k', 'LineWidth', 1.5)
    ylabel('Change in tremor severity')
    xlabel('Stimulation phase (degrees)')
    set(gca, 'XTickLabelRotation', 45)
    box('off')
    
    %raw ARC with Phase shift thereshold for main axes
%     f4= figure(4)
%     subplot(2,5,i)
%     bar(0:30:330,raw_main(i,:),'FaceColor',cl,'EdgeColor',cl)
%     hold on
%     yline(prctile(label_shift_1(i,:),99.7917),'k--','LineWidth',1)
%     yline(prctile(label_shift_1(i,:),0.2083),'k--','LineWidth',1)
%     ylabel('Change in tremor severity')
%     xlabel('Stimulation phase (degrees)')
%     set(gca,'XTickLabelRotation',45)
%     set(gca,'FontSize',12)
%     box('off')
    
end

f1.Units = 'centimeters'; f2.Units = 'centimeters'; f3.Units = 'centimeters'; f4.Units = 'centimeters';
f1.OuterPosition = [10, 10, 55, 15]; f2.OuterPosition = [10, 10, 50, 15]; f3.OuterPosition = [10, 10, 50, 15]; f4.OuterPosition = [10, 10, 50, 15];
set(f1, 'color', 'w'); set(f2, 'color', 'w'); set(f3, 'color','w');%set(f4,'color','w');

% %% Plot raw
% for i=1:length(cohort)
%     figure(4)
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

%raw ARC with NS thereshold for 3 axis

for i=1:size(tt1,1)
    f1 = figure()
    for axis = 1:size(tt1, 2)
        subplot(1, 3, axis)
        bar(0:30:330, squeeze(raw(i, axis,:)),'FaceColor',cl,'EdgeColor',cl)
        hold on
        yline(prctile(ns_mat1(i, axis, :), 99.7917),'k--','LineWidth',1)
        yline(prctile(ns_mat1(i, axis, :), 0.2083),'k--','LineWidth',1)
        
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



% % ARC median splot
% for hh=1:size(arc1,1)
%     f1=figure(1)
%     subplot(2,5,hh)
%     y=[arc1(hh,:);arc2(hh,:)]';
%     b = bar(0:30:330,y);
%     yline(0,'LineWidth',1)
% %     ylim([-0.75 0.75])
% %     yticks([ -1:0.25:1])
%     box('off')
%     ylabel({'Change in tremor severity'})
%     xlabel({'Stimulation phase (degrees)'})
%     
%     f1.Units = 'centimeters';
%     f1.OuterPosition= [10, 10, 50, 15];
%     set(gca,'XTickLabelRotation',45)
%     set(gca,'FontSize',12)
%     set(f1,'color','w');   
% end


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
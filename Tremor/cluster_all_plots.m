clear all
close all

c = 2; %% Cluster number - change it here


load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/clusters_BA.mat','nostim_c1')
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/clusters_BA.mat','nostim_c2')
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/output_auto.mat','TT1_C1')
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/output_auto.mat','TT1_C2')


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
y=[];
for i = 1:length(cohort)
   nostim_1(i,:) = squeeze(nostim(i, 1, :));
   if isempty(tt1{i,1}) || isnan(nostim_1(i, 1))
       idx = [idx, i];
   else
       y=[y i];
   end
end
clear i nostim_1

tt1=tt1(y,:);
cohort=cohort(y);
nostim=nostim(y,:,:);

clearvars -except c cl y tt1 cohort nostim

% change=[ z z x z x x x x z z];
% main=[1 1 3 1 3 3 3 3 1 1];
% ns_mat={[1 2 3];[1 2 3];[3 2 1];[1 2 3];[3 2 1];[3 2 1];[3 2 1];[3 2 1];[1 2 3];[1 2 3]};

%%%%%%%% INACABADO - MUDAR DE ACORDO COM OS RESULTADOS!
if c == 1
%     if strcmp(method, 'ward-pca2')     
        axi = [NaN 3 1 1 NaN 1 3 NaN 1 2]; 
        main=[NaN 3 3 1 NaN 3 1 NaN 1 2]; 
        main=main(y);
        ns_mat = [repmat(NaN,1,3); [3 2 1]; [3 2 1]; [1 2 3];repmat(NaN,1,3); [3 2 1]; [1 2 3]; repmat(NaN,1,3); [1 2 3]; [2 3 1]]; 
        ns_mat=ns_mat(y,:);
%     end
end

if c == 2 %% MUST ADD CASE TO REMOVE FAULTY BASELINE
%     if strcmp(method, 'ward-pca2')     
        axi = [1 1 3 1 1 NaN 1 1 1 1]; 
        main = [1 1 1 1 3 NaN 3 3 1 1];
        main=main(y);
        ns_mat = [[1 2 3]; [1 2 3]; [1 2 3]; [1 2 3]; [3 2 1]; repmat(NaN,1,3); [3 2 1]; [3 2 1]; [1 2 3]; [1 2 3]];
        ns_mat=ns_mat(y,:);
%     end
end
%% Plots

for ik = 1:size(tt1,1)
    for kk = 1:size(tt1,2)
        raw(ik,kk,:) = nanmedian(tt1{ik,kk});
        ns_mat1(ik,kk,:) = squeeze(nostim(ik, ns_mat(ik,kk),:));
    end   
end

%raw ARC with NS thereshold for 3 axis
load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','squash');
cl=blushred;

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
    filename=['CI_clust',num2str(c)','_pt',num2str(y(i)),'.png'];
    saveas(gcf,filename)
    filename=['CI_clust',num2str(c)','_pt',num2str(y(i)),'.fig'];
    saveas(gcf,filename)
end
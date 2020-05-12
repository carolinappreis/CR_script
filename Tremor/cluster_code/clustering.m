
function [clust]=clustering(all,iii,runs,clust,pc_comp);

clust.alldata=all;

Z = linkage(all,'ward');
c = cluster(Z, 'Maxclust', 2);

clust.C{iii,1}=c(1:runs(1));
clust.C{iii,2}=c(runs(1)+1:runs(1)+runs(2));
clust.C{iii,3}=c(runs(1)+runs(2)+1:end);
for i=1:3
    clust.count{iii,1}(:,i)=[numel(find(clust.C{iii,i}==1)); numel(find(clust.C{iii,i}==2))];
    clust.percent{iii,1}(:,i)=[numel(find(clust.C{iii,i}==1))./runs(i)*100; numel(find(clust.C{iii,i}==2))./runs(i)*100];
end

clust.win(iii,1)=find([clust.count{iii,1}(1,2) clust.count{iii,1}(2,2)]==(max([clust.count{iii,1}(1,2) clust.count{iii,1}(2,2)])));

figure(1)
subplot(2,5,iii)
[p h]=silhouette(all(runs(1)+1:runs(1)+runs(2),:),clust.C{iii,2});

for i=1:2
    clus=[];
    clus=find(clust.C{iii,2}==i);
    silh(1,i)=((numel(find((p(clus))<0)))./length(clus))*100;
    clust.mslh(iii,i)=mean(p(clus));
end

if ~isnan(clust.mslh(iii,clust.win(iii))) & clust.mslh(iii,clust.win(iii))<0.5 
    error('bad silhouette score for main cluster')
end

clust.idx{iii,1}=find(clust.C{iii,1}==clust.win(iii,1));
clust.idx{iii,2}=find(clust.C{iii,2}==clust.win(iii,1));
clust.idx{iii,3}=find(clust.C{iii,3}==clust.win(iii,1));

% figure(2);
% subplot(2,5,iii)
% bar(clust.percent{iii,1},'EdgeColor','none')
% box('off')
% legend({'NS','RS','PLS'},'Location','northwest')
% legend('boxoff')
% if clust.win(iii)==1
%     xticklabels({'cluster1*','cluster2'})
% else
%     xticklabels({'cluster1','cluster2*'})
%     ylabel('percentage of trials')
% end

% figure(3)
% set(gcf, 'color', 'w', 'Position', [300,300,1200,300])
% 
% subplot(1,3,1)
% scatter3([pc_trials_ns(:, 1); pc_trials(:, 1)], [pc_trials_ns(:, 2); pc_trials(:, 2)], [pc_trials_ns(:, 3); pc_trials(:, 3)], 10, c)
% title(sprintf('Patient %d - Combined', cohort(iii)))
% 
% subplot(1,3,2)
% scatter3(pc_trials_ns(1:5e4, 1), pc_trials_ns(1:5e4, 2), pc_trials_ns(1:5e4, 3), 10, c(1:5e4))
% title('Resampled Baseline')
% 
% subplot(1,3,3)
% scatter3(pc_trials(:, 1), pc_trials(:, 2), pc_trials(:, 3), 10, c_rs)
% title('Random Stim')

end
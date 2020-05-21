
function [clust]=clustering2(all,iii,runs,clust,pc_comp)


%%%----- cluster analysis
clust.alldata{iii,1}=all(1:runs(1)+runs(2),:);

Z = linkage(clust.alldata{iii,1},'ward');
c = cluster(Z, 'Maxclust', 2);

clust.C{iii,1}=c(1:runs(1));
clust.C{iii,2}=c(runs(1)+1:end);


for i=1:2
    clust.count{iii,1}(:,i)=[numel(find(clust.C{iii,i}==1)); numel(find(clust.C{iii,i}==2))];
    clust.percent{iii,1}(:,i)=[numel(find(clust.C{iii,i}==1))./runs(i)*100; numel(find(clust.C{iii,i}==2))./runs(i)*100];
end



%%% ------ silhouette analysis
figure(1)
subplot(2,5,iii)
[l h]=silhouette(all(1:runs(1),:),clust.C{iii,1});
p{1,1}=l;clear l h
figure(2)
subplot(2,5,iii)
[l h]=silhouette(all(runs(1)+1:runs(1)+runs(2),:),clust.C{iii,2});
p{2,1}=l;clear l h


for h=1:size(clust.C,2)
    for i=1:2
        clus=find(clust.C{iii,h}==i);
        clust.mslh(iii,h,i)=mean(p{h,1}(clus)); clear clus
    end
end

%%% ----- decision tree for which cluster vs non-cluster

mr=[clust.count{iii,1}(1,2) clust.count{iii,1}(2,2)];
if mr(1)~=mr(2)
    clust.win(iii,1)=find(mr==(max(mr)));
else
    pr=[clust.mslh(iii,2,1)  clust.mslh(iii,2,2)];
    clust.win(iii,1)=find(pr==(max(pr)));
end


if clust.mslh(iii,1,clust.win(iii,1))>0.75 
    clust.win(iii,1)=clust.win(iii,1);
    
    clust.idx{iii,1}=find(clust.C{iii,1}==clust.win(iii,1));
    clust.idx{iii,2}=find(clust.C{iii,2}==clust.win(iii,1));

else
    clust.win(iii,1)=NaN;
    
    clust.idx{iii,1}=1:length(clust.C{iii,1});
    clust.idx{iii,2}=1:length(clust.C{iii,2});
end



figure(4);
subplot(2,5,iii)
bar(clust.percent{iii,1},'EdgeColor','none')
box('off')
legend({'NS','RS'},'Location','northwest')
legend('boxoff')
if clust.win(iii)==1
    xticklabels({'cluster1*','cluster2'})
else
    xticklabels({'cluster1','cluster2*'})
    ylabel('percentage of trials')
end


% if ~isnan(clust.mslh(iii,clust.win(iii))) & clust.mslh(iii,clust.win(iii))<0.5
%     error('bad silhouette score for main cluster')
% end



figure(5)
set(gcf, 'color', 'w', 'Position', [300,300,1200,300])

subplot(1,4,1)
scatter3([pc_comp{iii,1}(:, 1); pc_comp{iii,2}(:, 1) ],[pc_comp{iii,1}(:, 2); pc_comp{iii,2}(:, 2) ],[pc_comp{iii,1}(:, 3); pc_comp{iii,2}(:, 3)],10,c);
title('Combined')

subplot(1,4,2)
scatter3(pc_comp{iii,1}(:, 1), pc_comp{iii,1}(:, 2), pc_comp{iii,1}(:, 3), 10, clust.C{iii,1})
title('Resampled Baseline')

subplot(1,4,3)
scatter3(pc_comp{iii,2}(:, 1), pc_comp{iii,2}(:, 2), pc_comp{iii,2}(:, 3), 10, clust.C{iii,2})
title('Random Stim')

end
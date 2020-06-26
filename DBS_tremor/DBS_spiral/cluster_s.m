
function [clust]= cluster_s(pc_t,iii)

x_all{iii,1}=[pc_t{iii,1}; pc_t{iii,2}];
runs{1,1}=1:length(pc_t{iii,1});
runs{1,2}=length(pc_t{iii,1})+1:length(x_all{iii,1});
k=3;
Z = linkage(x_all{iii,1},'ward');
c = cluster(Z, 'Maxclust', k);

figure()
subplot(3,1,1)
scatter3([pc_t{iii,1}(:, 1); pc_t{iii,2}(:, 1)], [pc_t{iii,1}(:, 2); pc_t{iii,2}(:, 2)], [pc_t{iii,1}(:, 3); pc_t{iii,2}(:, 3)], 10, c)
xlabel('x axis');ylabel('y axis');zlabel('z axis')
subplot(3,1,2)
scatter3([pc_t{iii,1}(:, 1)], [pc_t{iii,1}(:, 2)], [pc_t{iii,1}(:, 3)], 10, c(runs{1,1}))
xlabel('x axis');ylabel('y axis');zlabel('z axis')
subplot(3,1,3)
scatter3([pc_t{iii,2}(:, 1)], [pc_t{iii,2}(:, 2)], [pc_t{iii,2}(:, 3)], 10, c(runs{1,2}))
xlabel('x axis');ylabel('y axis');zlabel('z axis')
        
figure()
subplot(3,1,1)
[l h]=silhouette(x_all{iii,1},c);
title('ALL')
subplot(3,1,2)
[l h]=silhouette(pc_t{iii,1},c(runs{1,1}));
title('NS')
subplot(3,1,3)
[l h]=silhouette(pc_t{iii,2},c(runs{1,2}));
title('RS')


for i=1:2
    clust.C{iii,i}=c(runs{1,i});
    clust.idx{iii,i}=[{find(clust.C{iii,i}==1)}; {find(clust.C{iii,i}==2)}; {find(clust.C{iii,i}==3)}];
    clust.count{iii,1}(1:3,i)=[numel(find(clust.C{iii,i}==1)); numel(find(clust.C{iii,i}==2)); numel(find(clust.C{iii,i}==3))];
    clust.percent{iii,1}(1:3,i)=[numel(find(clust.C{iii,i}==1))./length(runs{1,i})*100; numel(find(clust.C{iii,i}==2))./length(runs{1,i})*100; numel(find(clust.C{iii,i}==3))./length(runs{1,i})*100];   
end

figure();
subplot(1,4,iii)
bar(clust.percent{iii,1},'EdgeColor','none')
box('off')
legend({'NS1','NS2','NS3'},'Location','northwest')
legend('boxoff')
ylabel('percentage of trials')



end


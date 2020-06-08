
function [clust,out]=clustering2(out,iii,clust,start,ending,yy)

load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/clustering_xc.mat','x_all')

%%%----- cluster analysis
all1=x_all{iii,1};
runs{1,1}=1:5e4;
runs{1,2}=5e4+1:size(all1,1);


Z = linkage(all1,'ward');
c = cluster(Z, 'Maxclust', 2);

for i=1:2
    clust.C{iii,i}=c(runs{1,i});
    clust.count{iii,1}(:,i)=[numel(find(clust.C{iii,i}==1)); numel(find(clust.C{iii,i}==2))];
    clust.percent{iii,1}(:,i)=[numel(find(clust.C{iii,i}==1))./length(runs{1,i})*100; numel(find(clust.C{iii,i}==2))./length(runs{1,i})*100];
    
    %%% ------ silhouette analysis
    figure(2+i)
    subplot(2,5,iii)
    [l h]=silhouette(all1(runs{1,i},:),clust.C{iii,i});
    p{i,1}=l;clear l h
end


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


load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/NS_all','bs_end','bs_begin')

co=1;
out.start_c{iii,co}=bs_begin(iii,clust.idx{iii,co});
out.ending_c{iii,co}=bs_end(iii,clust.idx{iii,co});
out.sns{iii,co}=start{iii,co};
out.ens{iii,co}=ending{iii,co};

co=2;
out.start_c{iii,co}=start{iii,co}(clust.idx{iii,co});
out.ending_c{iii,co}=ending{iii,co}(clust.idx{iii,co});
out.yy{iii,co}=yy{iii,co}(clust.idx{iii,co});

co=3;
out.start_c{iii,co}=start{iii,co};
out.ending_c{iii,co}=ending{iii,co};
out.yy{iii,co}=yy{iii,co};





% figure(5);
% subplot(2,5,iii)
% bar(clust.percent{iii,1},'EdgeColor','none')
% box('off')
% legend({'NS','RS'},'Location','northwest')
% legend('boxoff')
% ylabel('percentage of trials')
% if clust.win(iii)==1
%     xticklabels({'cluster1*','cluster2'})
% elseif clust.win(iii)==2
%     xticklabels({'cluster1','cluster2*'})
%
% else
%     xticklabels({'cluster1','cluster2'})
% end


% if ~isnan(clust.mslh(iii,clust.win(iii))) & clust.mslh(iii,clust.win(iii))<0.5
%     error('bad silhouette score for main cluster')
% end


% NOT FINISHED
% figure(6)
% set(gcf, 'color', 'w', 'Position', [300,300,1200,300])
%
% subplot(1,4,1)
% scatter3([pc_comp{iii,1}(:, 1); pc_comp{iii,2}(:, 1) ],[pc_comp{iii,1}(:, 2); pc_comp{iii,2}(:, 2) ],[pc_comp{iii,1}(:, 3); pc_comp{iii,2}(:, 3)],10,c);
% title('Combined')
%
% subplot(1,4,2)
% scatter3(pc_comp{iii,1}(:, 1), pc_comp{iii,1}(:, 2), pc_comp{iii,1}(:, 3), 10, clust.C{iii,1})
% title('Resampled Baseline')
%
% subplot(1,4,3)
% scatter3(pc_comp{iii,2}(:, 1), pc_comp{iii,2}(:, 2), pc_comp{iii,2}(:, 3), 10, clust.C{iii,2})
% title('Random Stim')

end

function [clust,out]=clustering2(out,iii,clust,start,ending,yy)

load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/aux_out.mat','x_all','pc1_exp')

%%%----- cluster analysis
all1=x_all{iii,1};
runs{1,1}=1:5e4;
runs{1,2}=5e4+1:size(all1,1);
k=2;

Z = linkage(all1,'ward');
c = cluster(Z, 'Maxclust', k);

for i=1:k
    clust.C{iii,i}=c(runs{1,i});
    clust.count{iii,1}(:,i)=[numel(find(clust.C{iii,i}==1)); numel(find(clust.C{iii,i}==2))];
    clust.percent{iii,1}(:,i)=[numel(find(clust.C{iii,i}==1))./length(runs{1,i})*100; numel(find(clust.C{iii,i}==2))./length(runs{1,i})*100];
    
    %%% ------ silhouette analysis
    figure(2+i)
    subplot(2,5,iii)
    [l h]=silhouette(all1(runs{1,i},:),clust.C{iii,i});
    xlim([-1 1])
%     xticklabels([-1:0.2:1])
    p{i,1}=l;clear l h
    title(sprintf('patient %d',(iii)))
    box('off')
end
out.bsilh{iii,1}=p{1,1};

for h=1:size(clust.C,2) %% baseline vs. random stim
    for i=1:k  %% cluster 1 vs. cluster 2
        clus=find(clust.C{iii,h}==i);
        clust.mslh(iii,h,i)=mean(p{h,1}(clus)); clear clus
    end
end

%%% ----- decision for which cluster vs non-cluster

rm=[clust.count{iii,1}(1,2) clust.count{iii,1}(2,2)];
if rm(1)~=rm(2)
    clust.win(iii,1)=find(rm==(max(rm)));
else
    bs=[clust.count{iii,1}(1,1) clust.count{iii,1}(2,1)];
    clust.win(iii,1)=find(bs==(max(bs)));
end


if clust.mslh(iii,1,clust.win(iii,1))>0.75 && clust.mslh(iii,2,clust.win(iii,1))>0.75 
    clust.win(iii,1)=clust.win(iii,1);
    for b=1:2
    clust.idx{iii,b}=find(clust.C{iii,b}==clust.win(iii,1));
    end
else
    clust.win(iii,1)=NaN;
    
    for b=1:2
    clust.idx{iii,b}=1:length(clust.C{iii,b});
    end
end


load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/aux_out.mat','bs_end','bs_begin')

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



figure(5);
subplot(2,5,iii)
bar(clust.percent{iii,1},'EdgeColor','none')
box('off')
legend({'NS','RS'},'Location','northwest')
legend('boxoff')
ylabel('percentage of trials')
ylim([0 100])
title(sprintf('patient %d',(iii)))
if clust.win(iii)==1
    xticklabels({'cluster1*','cluster2'})
elseif clust.win(iii)==2
    xticklabels({'cluster1','cluster2*'})

else
    xticklabels({'cluster1','cluster2'})
end
    
end
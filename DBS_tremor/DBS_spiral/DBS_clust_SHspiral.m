clear
load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/P01_NS1_SH.mat')
bins=find(NS1(:,1)==-1);
for y=1:length(bins)
    if(y+1)<length(bins)
        data1{1,y}=NS1(bins(y):bins(y+1),:);
    end
end
data1{1,length(bins)-1}=NS1(bins(y):end,:);

for nr=4:6
    data=data1{1,nr};
    
    L = 150;
    for o=1:3
        v=data(:,o);
        n = floor(length(v) / L); % Number of segments
        idx_st = NaN(n, 1);
        idx_st(1, 1) = 1;
        idx_end = NaN(n, 1);
        idx_end(1, 1) = L;
        v_new = NaN(n, L);
        for i = 2:n
            idx_st(i) = idx_st(i-1) + L;
            idx_end(i) = idx_end(i-1) + L;
        end
        for i = 1:n
            v_new(i, :) = v( idx_st(i, 1) : idx_end(i, 1),1);
        end
        resamp{1,o}=v_new'; clear v_new
    end
    for j = 1:size(resamp{1,1},1)
        x = [resamp{1,1}(j,:); resamp{1,2}(j,:)];
        [pc, score, latent, tsquare, explained] = pca(x');
        pc_trials_ns(j, 1:2) = pc(1:2, 1);
        explained_ns(j, 1:2) = explained;
        clear x
    end
    %         figure(4)
    plot(pc_trials_ns(:,1), 'r.')
    hold on
    plot(pc_trials_ns(:,2), 'b.')
    pc_t=pc_trials_ns;
    x_all=pc_t;
    runs=1:length(x_all);
    
    k=2;
    Z = linkage(x_all,'ward');
    c = cluster(Z, 'Maxclust', k);
    scatter([pc_t(:, 1)], [pc_t(:, 2)], 10,c)
    figure()
    [l h]=silhouette(x_all,c);
    
    for i=1:k
        clust.C=c(runs);
        clust.idx{i,1}=[find(clust.C==i)];
        clust.count(i,1)=[numel(find(clust.C==i))];
        clust.percent(i,1)=[numel(find(clust.C==i))./length(runs)*100];
    end
    
    p=[[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250]; [0 0 0]];
    figure (10)
    plot(data(:,1),data(:,2),'Color',[0.5 0.5 0.5])
    hold on
    
    for cl=1:k
        id=[clust.idx{cl,1}];
        dum=resamp{1,2}(id,:);
        tdum=resamp{1,1}(id,:);
        for ii=1:size(dum,1)
            plot(tdum(ii,:),dum(ii,:),'.','Color',p(cl,:),'MarkerSize',10)
            hold on
        end
        clear dum t dum
        clear id
    end
    close all
end
clear all
%  load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/all_clust.mat')
 load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/clustering_xc.mat')


for ii =1:size(x_all,1)
    
    data= x_all{ii,1}(5e4+1:end, 1:3);
    clust=c_all{ii,1}(5e4+1:end);
    f1=figure(1)
    subplot(2,5,ii)
    [p h]=silhouette(data,clust);
    
    for i=1:2
        clus=[];
        clus=find(clust==i);
        silh(1,i)=((numel(find((p(clus))<0)))./length(clus))*100;
        m_s(ii,i)=mean(p(clus));
    end
    
    E1=evalclusters(data,'linkage','Silhouette','klist',[1:2]) 
    % E=evalclusters(X,c,'Silhouette')
    optk(ii,1)=E1.OptimalK;
    critval(ii,:)=E1.CriterionValues;
    per_neg(ii,:)=silh;
    
    clearvars -except ii x_all c_all f1 optk per_neg critval m_s
    
end


% dendrogram(Z)
% Y=pdist(X);
% D = cophenet(Z,Y);  %'DavisBouldain' 'Gap' 'Silhouette' 'CalinskiHarabasz'
% Davies-Bouldin Index evaluates intra-cluster similarity and inter-cluster differences
%The Silhouette Index measure the distance between each data point, the centroid of the cluster it was assigned to and the closest centroid belonging to another cluster


data=[bfilt_up{1,1} ]';
data1= (cell2mat(filt_up{1,1}))';

Z = linkage(data, 'ward');
c = cluster(Z, 'Maxclust', 2);
scatter3(data(:, 1), data(:, 2), data(:, 3), 10, c)

load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/aux_out.mat')
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cluster_out_mc.mat');

for iii=1:10
    for a=1:3
    rs{iii,1}=x_all{iii,1}(5e4+1:end,:);
    end
end


for iii=1:10
    for a=1:3
    rs1{iii,1}=rs{iii,1}(clust.idx{iii,2},:);
    end
end

for iii=1:10
    for p=1:size(rs1{iii,1},1)
    main{iii,1}(p)= find(rs1{iii,1}(p,:)==max(rs1{iii,1}(p,:)));
    suma(iii,1:2)=[((numel(find(main{iii,1}==1))./size(rs1{iii,1},1).*100))  ((numel(find(main{iii,1}==2)) + numel(find(main{iii,1}==3)))./size(rs1{iii,1},1)).*100 ];
    end
    
    subplot(2,5,iii)
    pie(suma(iii,:))
   
end
bar(median(suma))
hold on
plot(1:2,suma','MarkerSize',15)
xlim([0 3])




p1=figure;
boxplot(suma)
box('off')
set(p1,'color','w')
ylabel('trials(%)')
xticklabels({'main axis','other axes'})

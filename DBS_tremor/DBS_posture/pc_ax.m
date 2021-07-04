function [f1]=pc_ax(s,tp_s,tp_e,cr,start,ending,color_b1,match_ax,pt)

dum=s.filt{1,3};
dum2=s.env_acc{1,3};

r=0;
for t=1:size(tp_s,1)
    for j = 1:length(tp_s{t,1})
        r=r+1;
        if (~isnan(tp_s{1,1}(j)))
            for ii=1:3
                x(ii,:) = dum(ii,(tp_s{t,1}(j)-3000):tp_s{t,1}(j));
                y(ii,:) = mean(dum2(ii,(tp_s{t,1}(j)-3000):tp_s{t,1}(j)));
            end
            [pc, score, latent, tsquare, explained] = pca(x');
            plot_pc(:,r)=pc(1:3,1);
            acc(:,r)=y;
            [k,kk]=sort(pc(1:3,1),'descend');
            pc_before(t,j) = kk(1);
            dif_before(t,j) = k(1)-k(2);
        end
        clearvars -except t j pc_before dum tp_s tp_e start dif_before ending cr plot_pc r acc dum2 color_b1 match_ax pt
    end
end



% plot(plot_pc','.')
% box('off')
% xlim([0 size(plot_pc,2)+1])


if cr==1
ending{1,3}(1,2)=258848;
end

dr=ending{1,3};
for j=1:size(dr,2)
    if (~isnan(dr(j)))
        for ii=1:3
            x(ii,:) = dum(ii,(dr(j)-3000):dr(j));
        end
        [pc, score, latent, tsquare, explained] = pca(x');
        [k,kk]=sort(pc(1:3,1),'descend');
        pc_last(1,j) = kk(1);
        dif_last(1,j)=k(1)-k(2);
    end
    clearvars -except t j pc_before dum tp_s tp_e pc_last cr dif_last dif_before start ending dr color_b1 match_ax pt
end


for t=1:size(tp_s,1)
    for j = 1:length(tp_s{t,1})
        if (~isnan(tp_s{1,1}(j)))
            for ii=1:3
                x(ii,:) = dum(ii,(tp_e{t,1}(j)):tp_e{t,1}(j)+3000);
            end
            [pc, score, latent, tsquare, explained] = pca(x');
            [k,kk]=sort(pc(1:3,1),'descend');
            pc_after(t,j) = kk(1);
            dif_after(t,j) = k(1)-k(2);
        end
        clearvars -except t j pc_after pc_before dum tp_s tp_e pc_last dif_before dif_last dif_after start ending color_b1 match_ax pt
    end
end



dr=start{1,3};
for j=1:size(dr,2)
    if (~isnan(dr(j)))
        for ii=1:3
            x(ii,:) = dum(ii,(dr(j)-3000):dr(j));
        end
        [pc, score, latent, tsquare, explained] = pca(x');
        [k,kk]=sort(pc(1:3,1),'descend');
        pc_first(1,j) = kk(1);
        dif_first(1,j)=k(1)-k(2);
    end
        clearvars -except pc_first dif_first cr t j pc_after pc_before dum tp_s tp_e pc_last dif_before dif_last dif_after start ending dr color_b1 match_ax pt
end

map=[pc_first' pc_before(:,1) pc_after(:,1) pc_before(:,2) pc_after(:,2) pc_last'];

map1=[ dif_first' dif_before(:,1) dif_after(:,1) dif_before(:,2) dif_after(:,2) dif_last'];

if match_ax(2,pt,1)==3
    mp=zeros(3,6);
  for i=1:3
      mp(i,find(map(i,:)==1))=3;
      mp(i,find(map(i,:)==3))=1;
  end
else
    mp=map;
end
      
f1=figure(4)
for i=1:size(mp,1)
    subplot(1,3,i)
    b=bar(mp(i,2:end),0.3,'FaceAlpha',[0.75],'EdgeColor','none')
    b(1).FaceColor=color_b1(2,:);
    xticklabels({'Btap','Atap','Btap','Atap','Boff'})
    xtickangle(45)
    xlim([0.8 5.2])
    xticks(1:5)
    ylabel('tremor axis')
    box('off')
    ylim([0 3])
    yticks([1:3])
    yticklabels({'Z','Y','X'})
     set(gca,'FontSize',11);
end
f1.Units = 'centimeters';
f1.OuterPosition= [10, 10, 35, 8];
set(f1,'color','w');



% figure(2)
% for i=1:3
%     subplot(1,3,i)
%     bar(map1(i,:))
%     xticklabels({'bef stim on';'bef tap';'after tap';'bef tap';'after tap';'bef stim off'})
%     xtickangle(45)
%     ylabel('difference of pca contrib')
%     box('off')
% end
end

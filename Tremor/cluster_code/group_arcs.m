function [stats]=group_arcs(out,stats)

load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','stone','squash');
color_b1=[stone ; blushred; aegean];
m_ax=1;
low=zeros(10,12);
high=zeros(10,12);
m_arc_a=zeros(10,12);
m_arc_s=zeros(10,12);

for iii=1:size(out.start_c,1)
    ref(iii,:)=nanmedian(squeeze(out.change_c{iii,2}{m_ax,1}));
    
    ref_ns(1:1e6,1:12)=squeeze(out.ns_arc{iii,m_ax});
    [low,high,m_arc_a,m_arc_s]=ns_th(out,ref_ns,iii,low,high,m_arc_a,m_arc_s); clear ref_ns
end
%%%% saved output from above
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/low_th_group')
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/high_th_group')

id_s=[];
id_a=[];
for i =1:size(ref,1)
    if (numel(find(ref(i,:)<0)))~=0 % all -except all amplifying subjects
        id_s=[id_s i];
    end
    if(numel(find(ref(i,:)>0)))~=0 % all -except all supressive subjects
        id_a=[id_a i];
    end
end

for cr=1:2
    if cr==1
        new_ref=ref(id_a,:);  %  amplifying effect;
        [p,q]=(max(new_ref'));
        phase_peak=q;
    else
        new_ref=ref(id_s,:);  %  supressive effect;
        [p,q]=(min(new_ref'));
        phase_peak=q;
    end
    
    for ii=1:size(new_ref,1)
        
        % main peak at 0 deg
        if phase_peak(ii)==1
            s_al(ii,:)=new_ref(ii,:);
            
        else
            s_al(ii,:)=[new_ref(ii,phase_peak(ii):end) new_ref(ii,1:phase_peak(ii)-1)];
            
            check(ii,1)=sum(diff([phase_peak(ii):size(s_al,2) 1:phase_peak(ii)-1]));
            
        end
        arc{cr,1}=s_al;
    end
    clear s_al
end


low1=low(id_s,:);

% for iii=1:10
% subplot(2,5,iii)
% bar(arc{1,1}(iii,:))
% end
% 
% figure
% subplot(2,1,1)
% plot(high','b')
% hold on
% plot(median(high),'b','lineWidth',4)
% plot(arc{1,1}','r')
% plot(median(arc{1,1}),'r','lineWidth',4)
% subplot(2,1,2)
% plot(low1','b')
% hold on
% plot(median(low1),'b','lineWidth',4)
% plot(arc{2,1}','r')
% plot(median(arc{2,1}),'r','lineWidth',4)


% for g=1:12
% r(1,g)=kstest(arc{1,1}(:,g)-high(:,g))==1;
% end
% sum(r)

for g=1:12
    [j,h]=signrank(arc{1,1}(:,g),(high(:,g)));
    stats.group_amp(g,2)=j;
    stats.group_amp(g,1)=h;clear j h
%     stats.group_amp(g,1)=h<(0.05/12); clear j h
    [j,h]=signrank(arc{2,1}(:,g),low1(:,g));
    stats.group_sup(g,2)=j;
    stats.group_sup(g,1)=h;clear j h
%     stats.group_sup(g,1)=h<(0.05/12); clear j h
    dif_a(:,g)=arc{1,1}(:,g)-high(:,g);
    dif_s(:,g)=arc{2,1}(:,g)-low1(:,g);
end

load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','stone','squash','sapphire','azure');
pp=[blushred; stone];

f1=figure(50)
subplot(1,2,1)
b=bar([median(arc{1,1}(:,1))  median(high(:,1))],'FaceColor','flat','FaceAlpha',.7,'EdgeColor','none','BarWidth',1);
b.CData(1,:) = pp(1,:);
b.CData(2,:) = pp(2,:);
hold on
plot([arc{1,1}(:,1)' ;  high(:,1)'],'k.')
ylim([-1 1])
xticklabels({'stim','non-stim'})
ylabel({'maximum amplification' ; 'of tremor'})
box('off')
set(gca,'FontSize',12);

subplot(1,2,2)
b=bar([median(arc{2,1}(:,1)) median(low1(:,1))],'FaceColor','flat','FaceAlpha',.7,'EdgeColor','none','BarWidth',1);
b.CData(1,:) = pp(1,:);
b.CData(2,:) = pp(2,:);
hold on
plot([arc{2,1}(:,1)' ;  low1(:,1)'],'k.')
ylim([-1 1])
xticklabels({'stim','non-stim'})
ylabel({'maximum suppression' ; 'of tremor'})
box('off')
set(f1,'color','w');
set(gca,'FontSize',12);

end


% % f1=figure(19)
% % for cr=1:2
% %     subplot(1,2,cr)
% %     % main peak at 180 deg
% %     dr=[arc{cr,1}(:,8:12) arc{cr,1}(:,1:7)];
% %     y2=nanmedian(dr,1);
% %     y1=prctile(dr,75);
% %     y3=prctile(dr,25);
% %     phase=1:12;
% %     phase1=[6:12 1:5];
% %     patch([phase fliplr(phase)], [y1 fliplr(y2)],[color_b1(1,:)],'FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
% %     hold on
% %     patch([phase fliplr(phase)], [y2 fliplr(y3)],[color_b1(1,:)],'FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
% %     % plot(y2,'.', 'MarkerSize',20,'Color',color_b1(p,:))
% %     stem(phase,y2,'.', 'LineWidth',4,'MarkerSize',20,'Color',color_b1(3,:))
% %     yline(0)
% %     
% %     if cr==1
% %         plot(phase1(find(stats.group_amp(:,1)==1)),y2(phase1(find(stats.group_amp(:,1)==1))),'.','MarkerSize',25,'Color',color_b1(2,:))
% %     elseif cr==2
% %         plot(phase1(find(stats.group_sup(:,1)==1)),y2(phase1(find(stats.group_sup(:,1)==1))),'.','MarkerSize',25,'Color',color_b1(2,:))
% %     end
% %     
% %     xlim([1 12])
% %     xticks([ 3 6 9 12])
% %     xticklabels({'-90','0','90','180'})
% %     ylabel({'Change in tremor severity'})
% %     if cr==1
% %         xlabel({'Alignment to phase with max';'tremor amplification'})
% %     else
% %         xlabel({'Alignment to phase with max';'tremor supression'})
% %     end
% %     ylim([-0.82 0.82])
% %     set(f1,'color','w');
% % end


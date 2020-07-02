function [stats]=group_aligned(out,stats)

load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','stone','squash');
color_b1=[stone ; blushred; aegean];
m_ax=1;

for iii=1:size(out.start_c,1)
    ref(iii,:)=nanmean(squeeze(out.change_c{iii,2}{m_ax,1}));
    ref_ns(iii,:)=squeeze(out.change_c{iii,1}(m_ax,:));
end

cr1=mean(prctile(ref_ns,99.7917,2));
cr2=mean(prctile(ref_ns, 0.2083,2));

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
        new_ref_ns=ref_ns(id_a,:);
        [p,q]=(max(new_ref'));
        phase_peak=q;
    else
        new_ref=ref(id_s,:);  %  supressive effect;
        new_ref_ns=ref_ns(id_s,:);
        [p,q]=(min(new_ref'));
        phase_peak=q;
    end
    
    for ii=1:size(new_ref,1)
        
        % main peak at 0 deg
        if phase_peak(ii)==1;
            s_al(ii,:)=new_ref(ii,:);
            
        else
            s_al(ii,:)=[new_ref(ii,phase_peak(ii):end) new_ref(ii,1:phase_peak(ii)-1)];
            
            check(ii,:)=[phase_peak(ii):size(s_al,2) 1:phase_peak(ii)-1];
            
        end
        arc{cr,1}=s_al;
    end
    f1=figure(19)
    subplot(2,1,cr)
    % main peak at 180 deg
    dr=[s_al(:,8:12) s_al(:,1:7)];
    y2=nanmedian(dr,1);
    y1=prctile(dr,75);
    y3=prctile(dr,25);
    time=1:12;
    patch([time fliplr(time)], [y1 fliplr(y2)],[color_b1(1,:)],'FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
    hold on
    patch([time fliplr(time)], [y2 fliplr(y3)],[color_b1(1,:)],'FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
    % plot(y2,'.', 'MarkerSize',20,'Color',color_b1(p,:))
    stem(time,y2,'.', 'LineWidth',4,'MarkerSize',20,'Color',color_b1(3,:))
    yline(0)
    h=find(y2>cr1 | y2<cr2);
    if ~isempty(h)
        plot(time(h),y2(h),'.','MarkerSize',25,'Color',color_b1(2,:))
    end
%             yline(cr2,'--')
%             yline(cr1,'--')
    xlim([1 12])
    xticks([ 3 6 9 12])
    xticklabels({'-90','0','90','180'})
    ylabel({'Change in tremor severity'})
    if cr==1
        xlabel({'Alignment to phase with max';'tremor amplification'})
    else
        xlabel({'Alignment to phase with max';'tremor supression'})
    end
    ylim([-0.82 0.82])
    set(f1,'color','w');

    
    %     clearvars -except ref ref_ns id_a id_s res color_b1 f1 cr ref_nostim cr2 cr1 arc
end

% kstest(arc1-dr1)==1

dr1=prctile(ref_ns(id_a,:),95,2);
dr2=prctile(ref_ns(id_s,:), 5,2);

for g=1:size(arc{1,1},2)
    arc1=arc{1,1}(id_a,g);
    arc2=arc{2,1}(id_s,g);
    [j,h]=signrank(arc1,dr1);
    stats.group_amp(1,g)=j; 
    stats.group_amp(2,g)=j<(0.05/12);clear j h
    [j,h]=signrank(arc2,dr2);
    stats.group_sup(1,g)=j;
    stats.group_sup(2,g)=j<(0.05/12);
    clear arc1 arc2 j h
end

end


clear all; close all;
load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','squash','aegean','stone');
color_b1=[stone ; blushred];

%%% non cluster data!
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cleaned_rc12_noaddon.mat')
% load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/newnonstim10.mat')
% main=[1 1 3 1 3 3 3 3 1 1];
% for i=1:10
%     ref(i,:)=nanmedian(tt1{i,1});
%     ref_ns=squeeze((nostimout(:,main(i),:)));
% end



%%% clustered data
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/data_matchcluster')
ref=squeeze((tt1(:,1,:)));
ref_ns=squeeze((nostimout(:,1,:)));


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
            ns_al(ii,:)=new_ref_ns(ii,:);
            
        else
            s_al(ii,:)=[new_ref(ii,phase_peak(ii):end) new_ref(ii,1:phase_peak(ii)-1)];
            ns_al(ii,:)=[new_ref_ns(ii,phase_peak(ii):end) new_ref_ns(ii,1:phase_peak(ii)-1)];
            
            check(ii,:)=[phase_peak(ii):size(s_al,2) 1:phase_peak(ii)-1];
            
        end
        
        
        %%% 180o from main peak at 0 deg
        if (phase_peak(ii)+5)<= size(s_al,2)
            
            a_s2(ii,:)= [new_ref(ii,phase_peak(ii)+5:end) new_ref(ii,1:phase_peak(ii)+5-1)];
            a_ns2(ii,:)= [new_ref_ns(ii,phase_peak(ii)+5:end) new_ref_ns(ii,1:phase_peak(ii)+5-1)];
        else
            
            dum=(phase_peak(ii)+5)-size(s_al,2);
            a_s2(ii,:)= [new_ref(ii,dum+5:end) new_ref(ii,1:dum+5-1)];
            a_ns2(ii,:)= [new_ref_ns(ii,dum+5:end) new_ref_ns(ii,1:dum+5-1)];
            clear dum
            
        end
        
        
    end
    
    %     if kstest(a_s_al(:,1)-a_ns_al(:,1))==1 % wilcoxon else ttest
    %     because amp requires wilcoxon we left it as wilcoxon for both
    for i=1:12
        [p1,h1]=signrank(s_al(:,i),ns_al(:,i));
        ref_nostim(cr,i)=p1;
        
        [p2,h2]=signrank(s_al(:,i),a_s2(:,i));
        ref_180(1,i)=p2;
        
        [p3,h3]=signrank(ns_al(:,i),a_ns2(:,i));
        ref_ns_180(1,i)=p3;
    end
    
    cond={ns_al;s_al};
    
    f1=figure(1)
    
    subplot(1,2,cr)
    for m=1:2
        
        rr=cond{m,1};
        % main peak at 180 deg
        dr=[rr(:,8:12) rr(:,1:7)];
        phase=[ref_nostim(cr,8:12) ref_nostim(cr,1:7)];
        y2=nanmedian(dr,1);
        y1=prctile(dr,75);
        y3=prctile(dr,25);
        time=1:12;
        patch([time fliplr(time)], [y1 fliplr(y2)],[color_b1(m,:)],'FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
        hold on
        patch([time fliplr(time)], [y2 fliplr(y3)],[color_b1(m,:)],'FaceAlpha',[0.15],'EdgeColor','none','HandleVisibility','off')
        % plot(y2,'.', 'MarkerSize',20,'Color',color_b1(p,:))
        stem(time,y2,'.', 'LineWidth',1.5,'MarkerSize',12,'Color',color_b1(m,:))
        yline(0)
        if ~isempty(find(phase<(0.05/12)))
            %     plot(time(find(phase<(0.05/12))),y2(find(phase<(0.05/12))),'.','Color',color_b1(1,:),'MarkerSize',10)
            plot(time(find(phase<(0.05/12))),(-0.45),'k*','MarkerSize',5)
        end
        xlim([1 12])
        xticks([ 3 6 9 12])
        ylim([-0.5 0.5])
        yticks([-0.5:0.25:0.5])
        xticklabels({'-90','0','90','180'})
        ylabel({'Change in tremor severity'})
        if cr==1
            xlabel({'Alignment to phase with max';'tremor amplification'})
        else
            xlabel({'Alignment to phase with max';'tremor supression'})
        end
        set(gca,'FontSize',12)
    end
    clearvars -except ref ref_ns id_a id_s res color_b1 f1 cr ref_nostim
    
end


% % %old way of plotting
% %
% %   if cr==1
% %         cl=blushred;
% %         cl1=squash;
% %
% %     else
% %         cl=aegean;
% %         cl1=stone;
% %     end
% %
% %     f1=figure(1)
% %     subplot(1,2,cr)
% %     plot(metric1','Color',cl1)
% %     ylim([-1 1])
% %     xlim([1 12])
% %     xticks([ 3 6 9 12])
% %     xticklabels({'-90','0','90','180'})
% %     box('off')
% %     hold on
% %     plot(mean(metric1),'k','LineWidth',4,'Color', cl)
% %     yline(0,'k', 'LineWidth',1,'LineStyle','--')
% %     ylabel({'Change in tremor severity'})
% %     xlabel({'Alignment to phase with max';'tremor supression'})
% %     set(gca,'FontSize',14)
% %
% %

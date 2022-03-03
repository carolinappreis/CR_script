
function[total]=coef_share(s,start,ending,tp_s,tp_e,match_ax,pt)

%%%% you need 1 second (1000 ponits segments)


dum=s.filt{1,3};
dum2=s.env_acc{1,3};
[ref]=pt_post(dum2, start, ending,tp_s,tp_e);


for t=1:size(ref,1)
    r=0;
    for j = 2:2:length(ref{t,1})
        r=r+1;
        if ~isempty(length(ref{t,1}))
            tempo=floor((numel(ref{t,1}(j-1):ref{t,1}(j))/1000))*1000;
            period=ref{t,1}(j)-tempo:ref{t,1}(j);
            bins=1:1000:length(period);
            for g=1:length(bins)
                if (g+1)<length(bins)
                for ax=1:size(dum,1)
                    x(ax,g,:) = dum(ax,period(bins(g):bins(g+1)));
                    y(ax,g) = sum(dum2(ax,period(bins(g):bins(g+1))));
                end
                end
            end
        end
      
        for hh = 1:size(x,2)
            % in pca, rows are observations and columns are variables
            for_pca = squeeze(x(:,hh,:)); % should be 3 vs length(segments)
            [pc] = pca(for_pca');
            gha=pc(1:3, 1);
            accum(hh,:)=find(gha==(max(gha)));
            clear for_pca pc gha
        end
        total{1,t}(r,:)=[numel(find(accum==1))./length(accum) numel(find(accum==2))./length(accum) numel(find(accum==3))./length(accum)];
        clearvars -except s samplerate ref color_b1 match_ax total t j dum dum2 r pt
    end
    
    if match_ax(2,pt,1)~=3
        share=total{1,t}*100;
    else
        share=[total{1,t}(:,3) total{1,t}(:,2) total{1,t}(:,1)]*100;
    end
    
%     f1=figure(1)
%     subplot(1,3,t)
%     bar(1:size(total{1,t},1),share,0.4,'stacked')
%     xticklabels({'Btap','Atap','Btap','Atap','Boff'})
%     xtickangle(45)
%     xlim([0.8 5.2])
%     xticks(1:5)
%     ylim([0 100])
%     yticks(0:25:100)
%     ylabel('Share of axes')
%     box('off')
%     f1.Units = 'centimeters';
%     f1.OuterPosition= [10, 10, 35, 8];
%     set(f1,'color','w')


    f1=figure(t)
    for u=1:3
    subplot(1,3,u)
    pie(share(u,:))
    end
    box('off')
    f1.Units = 'centimeters';
    f1.OuterPosition= [10, 10, 35, 8];
    set(f1,'color','w')

end
share= [mean([total{1,1} ; total{1,2}; total{1,3}]) ; std([total{1,1} ; total{1,2}; total{1,3}])]
end



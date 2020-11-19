function[fig]=sub_to_ctx(coh_filts,name)

for pr=1:size(coh_filts,1)
    ctx=coh_filts{pr,1}(1,:);
    thal=coh_filts{pr,1}(2:end,:);
    ctx_env=abs(hilbert(ctx));
    
    for ct=1:size(thal,1)
        ref=thal(ct,:);
        env=abs(hilbert(ref));
        env=smoothdata(env,'movmean',50);
        onset1=bursts(env);
        
        for o=1:size(onset1,1)
            onset=onset1{o,1};
            el=500;
            for jj=1:length(onset)
                if onset(jj)>el && onset(jj)+el<length(ctx_env)
                    al_seg(jj,:)=((ctx_env(onset(jj)-el:onset(jj)+el)-median(ctx_env(onset(jj)-el:onset(jj))))./median(ctx_env(onset(jj)-el:onset(jj))));
                else
                    al_seg(jj,1:1001)=NaN;
                    
                end
            end
            
            for j=1:length(ctx_env)/100
                idx_sur=randi([501,(length(ctx_env)-el)],1,1);
                surr_seg(j,:)= ((ctx_env(idx_sur-el:idx_sur+el)-median(ctx_env(idx_sur-el:idx_sur)))./median(ctx_env(idx_sur-el:idx_sur)));
            end
            
            burst_amp(pr,o,ct,:)=nanmedian(al_seg,1);
            surr_amp(pr,o,ct,:)=nanmedian(surr_seg,1);
            
            clear al_seg surr_seg onset
        end
        clear ref env onset1 
    end
    clear ctx thal ctx_env
end

clearvars -except burst_amp surr_amp el name

cond={'burst_amp' 'surr_amp'};
for co=1:size(cond,1)
    for co2=1:size(cond,2)
        clear dat
        dat1=eval(cond{co,co2});
        
        for  bst=1:2
            dum=squeeze(dat1(:,bst,:,300:800)); %%% CHOOSING [-200 300] TO PLOT
            join=[];
            for i=1:size(dum,1)
                dum1=squeeze(dum(i,:,:));
                dum2=dum1(any(dum1,2),:);
                join=[join ; dum2]; clear dum1 dum2
            end
            data{co,co2}(bst,:,:)=join; clear join
        end
    end
end


load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','squash','blood','sky','aegean');
if name=='bz'
    color_b={squash blood};
else
    color_b={sky aegean};
end

time=1:el+1;
site=([[15+1 15+1 15+2 15+2];[15+4 15+4 15+5 15+5]]);

for co=1:size(cond,1)
    fig=figure(co);
    for on=1:2
        st=NaN(1,length(time));
        clear A; A=squeeze(data{co,1}(on,:,:));
        clear B; B=squeeze(data{co,2}(on,:,:));
        hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
        beg=[];
        beg=find(st(1,:)<0.05 & st2(1,:)~=0);
        sig_rise_all=[];
        if ~isempty(beg)
            sig_rise_all=[beg(1) beg(end)];
        end
        clear st st2
        
        y2=100.*(mean(A));
        y1=100.*(mean(A)+std(A)./sqrt(size(A,1)));
        y3=100.*(mean(A)-std(A)./sqrt(size(A,1)));
%         y1= 100.*(mean(A)+std(A));
%         y3= 100.*(mean(A)-std(A));
        plot(time, y2,'LineStyle','-', 'LineWidth',1.5,'Color',[color_b{co,on}])
        hold on
        patch([time fliplr(time)], [y1 fliplr(y2)],[color_b{co,on}],'FaceAlpha',[0.2],'EdgeColor','none')
        patch([time fliplr(time)], [y2 fliplr(y3)],[color_b{co,on}],'FaceAlpha',[0.2],'EdgeColor','none')
        if ~isempty (beg)
            patch([sig_rise_all(1) sig_rise_all(2) sig_rise_all(2) sig_rise_all(1)],site(on,:),color_b{co,on},'EdgeColor','none')
            clear beg
        end
    %         xline(200,'--',{'burst onset'},'LabelOrientation','horizontal','Color',[0.5 0.5 0.5],'LineWidth',2)
        xline(200,'--','Color',[0.5 0.5 0.5],'LineWidth',2)
        box('off')
        xlim ([0 500])
        ylim ([-20 20])
        xticks([0:100:500])
        xticklabels ({'-200','-100','0','100','200','300'})
        fig.Units = 'centimeters';
        fig.OuterPosition= [10, 10, 8, 10];
        fig.Color='w';
        set(gca,'FontSize',12)
        ylabel('change in beta amplitude (%)');
        xlabel('Time (ms)');
        
        hold on
    end
end

end
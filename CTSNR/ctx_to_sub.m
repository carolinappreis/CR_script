function [fig]=ctx_to_sub(coh_filts,name)

for pr=1:size(coh_filts,1)
    ctx=coh_filts{pr,1}(1,:);
    thal=coh_filts{pr,1}(2:end,:);
    
    env=abs(hilbert(ctx));
%     env=smoothdata(env,'movmean',50);
%     [onset1,offset1]=bursts(env);
    
    [onset,offset]=bursts_aligned(env,ctx);
    onset1=onset; clear onset
    offset1=offset; clear offset
    
    for o=1:size(onset1,1)
        onset=onset1{o,1};
        el=500;
        for jj=1:length(onset)
            if onset(jj)>el && onset(jj)+el<length(env)
                cr1(jj,:)=env(onset(jj)-el:onset(jj)+el);
                ref_seg(jj,:)=((env(onset(jj)-el:onset(jj)+el)-median(env(onset(jj)-el:onset(jj))))./median(env(onset(jj)-el:onset(jj))));
            else
                ref_seg(jj,1:1001)=NaN;
            end
        end
        for ct=1:size(thal,1)
            thal_env=abs(hilbert(thal(ct,:)));
            thal_env=smoothdata(thal_env,'movmean',50);
            for jj=1:length(onset)
                if onset(jj)>el && onset(jj)+el<length(thal_env)
                    cr(jj,:)=thal_env(onset(jj)-el:onset(jj)+el);
                    al_seg(jj,:)=((thal_env(onset(jj)-el:onset(jj)+el)-median(thal_env(onset(jj)-el:onset(jj))))./median(thal_env(onset(jj)-el:onset(jj))));
                else
                    al_seg(jj,1:1001)=NaN;
                end
            end
            
            for j=1:length(env)/1000
                idx_sur=randi([501,(length(env)-el)],1,1);
                surr_seg(j,:)= ((thal_env(idx_sur-el:idx_sur+el)-median(thal_env(idx_sur-el:idx_sur)))./median(thal_env(idx_sur-el:idx_sur)));
            end
            
           burst_amp(pr,o,ct,:)=nanmedian(al_seg,1);
            surr_amp(pr,o,ct,:)=nanmedian(surr_seg,1);
            
            clear thal_env al_seg sur_seg
        end
        ref_amp(o,pr,:)=nanmedian(ref_seg,1);
        clear onset
    end
    clear ctx thal env onset1 ref_seg
end

clearvars -except burst_amp surr_amp ref_amp el name coh_filts

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
color_ref={[0.5 0.5 0.5]  [0 0 0]};
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
            patch([sig_rise_all(1) sig_rise_all(2) sig_rise_all(2) sig_rise_all(1)],site(on,:),color_ref{co,on},'EdgeColor','none')
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


fig=figure;
for on=1:2
    clear A;A=squeeze(ref_amp(on,:,300:800));
    y2=100.*(mean(A));
    y1=100.*(mean(A)+std(A)./sqrt(size(A,1)));
    y3=100.*(mean(A)-std(A)./sqrt(size(A,1)));
    %         y1= 100.*(mean(A)+std(A));
    %         y3= 100.*(mean(A)-std(A));
    plot(time, y2,'LineStyle','-', 'LineWidth',1.5,'Color',[color_ref{1,on}])
    hold on
    patch([time fliplr(time)], [y1 fliplr(y2)],[color_ref{1,on}],'FaceAlpha',[0.2],'EdgeColor','none')
    patch([time fliplr(time)], [y2 fliplr(y3)],[color_ref{1,on}],'FaceAlpha',[0.2],'EdgeColor','none')
end
%         xline(200,'--',{'burst onset'},'LabelOrientation','horizontal','Color',[0.5 0.5 0.5],'LineWidth',2)
xline(200,'--','Color',[0.5 0.5 0.5],'LineWidth',2)
box('off')
xlim ([0 500])
if max(mean(A).*100)>50
    ylim ([-150 150])
else
    ylim ([-20 20])
end
xticks([0:100:500])
xticklabels ({'-200','-100','0','100','200','300'})
fig.Units = 'centimeters';
fig.OuterPosition= [10, 10, 8, 10];
fig.Color='w';
set(gca,'FontSize',12)
ylabel('change in beta amplitude (%)');
xlabel('Time (ms)');


% st=NaN(1,length(time));
% clear A; A=squeeze(data{co,1}(1,:,:));
% clear B; B=squeeze(data{co,1}(2,:,:));
% hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
% beg=[];
% beg=find(st(1,:)<0.05 & st2(1,:)~=0);
% sig_rise_all=[];
% if ~isempty(beg)
%     sig_rise_all=[beg(1) beg(end)];
% end
% clear st st2
% fig=figure
% for on=1:2
%     clear A; A=squeeze(data{1,1}(on,:,:));
%     y2=100.*(mean(A));
%     y1=100.*(mean(A)+std(A)./sqrt(size(A,1)));
%     y3=100.*(mean(A)-std(A)./sqrt(size(A,1)));
%     %         y1= 100.*(mean(A)+std(A));
%     %         y3= 100.*(mean(A)-std(A));
%     plot(time, y2,'LineStyle','-', 'LineWidth',1.5,'Color',[color_b{co,on}])
%     hold on
%     patch([time fliplr(time)], [y1 fliplr(y2)],[color_b{co,on}],'FaceAlpha',[0.2],'EdgeColor','none')
%     patch([time fliplr(time)], [y2 fliplr(y3)],[color_b{co,on}],'FaceAlpha',[0.2],'EdgeColor','none')
%     %         xline(200,'--',{'burst onset'},'LabelOrientation','horizontal','Color',[0.5 0.5 0.5],'LineWidth',2)
%     xline(200,'--','Color',[0.5 0.5 0.5],'LineWidth',2)
%     box('off')
%     xlim ([0 500])
%     ylim ([-20 20])
%     xticks([0:100:500])
%     xticklabels ({'-200','-100','0','100','200','300'})
%     fig.Units = 'centimeters';
%     fig.OuterPosition= [10, 10, 8, 10];
%     fig.Color='w';
%     set(gca,'FontSize',12)
%     ylabel('change in beta amplitude (%)');
%     xlabel('Time (ms)');
%     
%     hold on
% end
%     if ~isempty (sig_rise_all)
%         patch([sig_rise_all(1) sig_rise_all(2) sig_rise_all(2) sig_rise_all(1)],site(on,:),'r','EdgeColor','none')
% %         clear beg
%     end




end

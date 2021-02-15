function [fig]=ctx_to_sub_psi(coh_filts,name)
r=0;
for pr=1:size(coh_filts,1)
    ctx=coh_filts{pr,1}(1,:);
    thal=coh_filts{pr,1}(2:end,:);
    
    env=abs(hilbert(ctx));
%         env=smoothdata(env,'movmean',50);
%         [onset1,offset1]=bursts(env);

        [onset,offset]=bursts_aligned(env,ctx);
        onset1=onset; clear onset
        offset1=offset; clear offset
    
    ctx_phase=angle(hilbert(ctx));
    
    for o=1:size(onset1,1)
        clear onset
        onset=onset1{o,1};
        
        for ct=1:size(thal,1)
            clear thal_phase non_norm
            thal_phase=angle(hilbert(thal(ct,:)));
            non_norm=wrapToPi(ctx_phase-thal_phase);
            
            clear epochs_idx epochs_t
            for jj=1:length(onset)
                el=500;
                if onset(jj)>el && onset(jj)+el<length(non_norm)
                    epochs_idx(jj,:)=onset(jj)-el:onset(jj)+el;
                    epochs_t(jj,:)=non_norm(onset(jj)-el:onset(jj)+el);
                end
                
            end
            
            epochs_idx = epochs_idx(any(epochs_idx,2),:);
            epochs_t = epochs_t(any(epochs_t,2),:);
            
            clear epochs_idx_sur epochs_t_sur
            for n=1:(length(non_norm)/1000)
                idx_sur=randi([el+1,(length(non_norm)-el)],1,1);
                epochs_idx_sur(n,:)= idx_sur-el:idx_sur+el;
                epochs_t_sur(n,:)= non_norm(idx_sur-el:idx_sur+el);
            end
            
            ol=50;
            for z= 1:size(epochs_idx,1)
                for w=1:size(epochs_idx,2)
                    if  epochs_idx(z,w)+ol<length(non_norm)
                        ep_b1(z,w)=circ_r((non_norm(epochs_idx(z,w):epochs_idx(z,w)+ol))');
                        ep_t1(1,w)=circ_r(epochs_t(:,w));
                    end
                end
            end
            
            for z= 1:size(epochs_idx_sur,1)
                for w=1:size(epochs_idx_sur,2)
                    if  epochs_idx_sur(z,w)+ol<length(non_norm)
                        ep_b1_s(z,w)=circ_r((non_norm(epochs_idx_sur(z,w):epochs_idx_sur(z,w)+ol))');
                        ep_t1_s(1,w)=circ_r(epochs_t_sur(:,w));
                    end
                end
            end
            ep_t(pr,o,ct,:)=zscore(ep_t1); 
            ep_t_s(pr,o,ct,:)=zscore(ep_t1_s); clear ep_t1_s
            ep_b(pr,o,ct,:)=zscore(nanmean(ep_b1,1)); 
            ep_b_s(pr,o,ct,:)=zscore(nanmean(ep_b1_s,1)); clear ep_b1_s
            r=r+1;
            base_t(r,:)=mean(ep_t1(300:500)); clear ep_t1
            base_b(r,:)=mean(ep_b1(300:500)); clear ep_b1
            
        end
    end
    clearvars -except ep_b ep_b_s ep_t ep_t_s pr coh_filts el name base_t base_b r
end


cond={ 'ep_b' 'ep_b_s' ; 'ep_t'  'ep_t_s'};
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


load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','squash','blood','sky','aegean','grey');
if name=='bz'
    color_b={squash blood};
else
    color_b={sky aegean};
end

time=1:el+1;
site=([[2+0.3 2+0.3 2+0.4 2+0.4];[2+0.6 2+0.6 2+0.7 2+0.7]]);

color_ref={[0.5 0.5 0.5]  [0 0 0]};

baseline=round([mean(base_b) mean(base_t)],2);


for co=1:size(cond,1)
    fig=figure(co);
    for on=1:2
        st=NaN(1,length(time));
        clear A; A=squeeze(data{co,1}(on,:,:));
        clear B; B=squeeze(data{co,2}(on,:,:));
        sig_rise_all=[];
        beg=[];
        hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
        beg=find(st(1,:)<0.05 & st2(1,:)~=0);
        if ~isempty(beg)
            sig_rise_all=[beg(1) beg(end)];
        end
        clear st st2
        
        y2=(mean(A));
        y1=(mean(A)+std(A)./sqrt(size(A,1)));
        y3=(mean(A)-std(A)./sqrt(size(A,1)));
        % y1= (mean(A)+std(A));
        % y3= (mean(A)-std(A));
        plot(time, y2,'LineStyle','-', 'LineWidth',1.5,'Color',[(color_ref{1,on}+color_b{1,on})./2])
        hold on
        patch([time fliplr(time)], [y1 fliplr(y2)],[(color_ref{1,on}+color_b{1,on})./2],'FaceAlpha',[0.2],'EdgeColor','none')
        patch([time fliplr(time)], [y2 fliplr(y3)],[(color_ref{1,on}+color_b{1,on})./2],'FaceAlpha',[0.2],'EdgeColor','none')
        if ~isempty (beg)
            patch([sig_rise_all(1) sig_rise_all(2) sig_rise_all(2) sig_rise_all(1)],site(on,:),color_ref{1,on},'EdgeColor','none')
        end
       %         xline(200,'--',{'burst onset'},'LabelOrientation','horizontal','Color',[0.5 0.5 0.5],'LineWidth',2)
        xline(200,'--','Color',[0.5 0.5 0.5],'LineWidth',2)
        box('off')
        xlim ([0 500])
        ylim ([-3 3])
        xticks([0:100:500])
        xticklabels ({'-200','-100','0','100','200','300'})
        fig.Units = 'centimeters';
        fig.OuterPosition= [10, 10, 8, 10];
        fig.Color='w';
        set(gca,'FontSize',12)
        if co==1
        ylabel('PSI within bursts (zscore)');
        else
         ylabel('PSI across bursts (zscore)');
        end
        xlabel('Time (ms)');
        hold on
    end
        str = ['psi:',num2str(baseline(co))];
        t = text(-1.5,-2.5,str);
        s = t.FontSize;
        t.FontSize = 12;
end

end

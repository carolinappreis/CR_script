clear all
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/final_mats')
load('SNR_bua.mat');

%filtering probe signals and cortex in the beta band (+-5Hz peak coherence) if the summed coherence
%between 15-35Hz is more than 10% of the coherence in all freqs.
m=1;
samprate=1000;
for pr=1:size(bua,1)
    for ct=2:size(bua{pr,1},1)
        [Pxx_ind,F_ind]=mscohere(bua{pr,1}(1,:),bua{pr,1}(ct,:),samprate,[],samprate,samprate);
        Pxx_ind_beta=Pxx_ind(16:36);
        coher(ct-1,:)=Pxx_ind_beta;
        coh_cri(ct-1,1)=(sum(Pxx_ind_beta)/sum(Pxx_ind(1:end)));
        clear Pxx_ind F_ind
    end
    
    fool=max(coher');
    ref=find(fool==max(fool));
    if coh_cri(ref)>0.1
        filtrange=14+find(Pxx_ind_beta(ref,:)==max(Pxx_ind_beta(ref,:)));
        [b,a]=butter(2,[(filtrange-5)/(0.5*samprate) (filtrange+5)/(0.5*samprate)],'bandpass');
        thal_ctx(m,:,:)=[filtfilt(b,a,bua{pr,1}(ref+1,:)) ; filtfilt(b,a,bua{pr,1}(1,:))];
        m=m+1;
    end
    
    clear coher coh_cri
end

% ctx=squeeze(thal_ctx(:,1,:));
% [env_var]=burst_var(ctx)


for pr=1:size(thal_ctx,1)
    
    env_ctx=squeeze(abs(hilbert(thal_ctx(pr,2,:))));
    env_ctx=smoothdata(env_ctx,'movmean',50)';
    
    env=squeeze(abs(hilbert(thal_ctx(pr,1,:))));
    env=smoothdata(env,'movmean',50)';
    
    onset1=bursts(env);
    
    for o=1:size(onset1,1)
        onset=onset1{o,1};
        for jj=1:length(onset)
            if onset(jj)>500 && onset(jj)+500<length(env_ctx)
                al_seg(jj,:)=((env_ctx(onset(jj)-500:onset(jj)+500)-median(env_ctx(onset(jj)-500:onset(jj))))./median(env_ctx(onset(jj)-500:onset(jj))));
            end
        end
        burst_amp(o,pr,:)=nanmedian(al_seg);
    end
    
    for j=1:1000;
        idx_sur=randi([501,(length(env)-500)],1,1);
        surr_seg(j,:)= ((env_ctx(idx_sur-500:idx_sur+500)-median(env_ctx(idx_sur-500:idx_sur)))./median(env_ctx(idx_sur-500:idx_sur)));
    end
    surr_amp(pr,:)=nanmedian(surr_seg);
    
    
    clear env_ctx env_thal surr_seg al_seg
end


clearvars -except burst_amp surr_amp

time=1:1001;
color_b(1,:)= [0 0 0.5];
color_b(2,:)= [0.5 0 0.5];
site(1,:)=[14+0.1 14+0.1 14+0.9 14+0.9];
site(2,:)=[16+0.1 16+0.1 16+0.9 16+0.9];


fig=figure()
for on=1:2;
    st=NaN(1,1001);
    clear A; A=squeeze(burst_amp(on,:,:));
    clear B; B=surr_amp;
    hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
    beg=find(st(1,:)<0.05 & st2(1,:)~=0);
    if ~isempty(beg)
        sig_rise_all=[beg(1) beg(end)];
    end
    clear st st2
    
    y2=100.*(mean(A));
    y1=100.*(mean(A)+std(A)./sqrt(size(A,1)));
    y3=100.*(mean(A)-std(A)./sqrt(size(A,1)));
    % y1= 100.*(mean(A)+std(A));
    % y3= 100.*(mean(A)-std(A));
    p2=plot(time, y2,'LineStyle','-', 'LineWidth',1.5,'Color',color_b(on,:))
    patch([time fliplr(time)], [y1 fliplr(y2)],[color_b(on,:)],'FaceAlpha',[0.2],'EdgeColor','none')
    patch([time fliplr(time)], [y2 fliplr(y3)],[color_b(on,:)],'FaceAlpha',[0.2],'EdgeColor','none')
    if ~isempty (beg)
        patch([sig_rise_all(1) sig_rise_all(2) sig_rise_all(2) sig_rise_all(1)],site(on,:),color_b(on,:),'EdgeColor','none')
        clear beg
    end
    xline(500,'--',{'burst onset'},'LabelOrientation','horizontal','Color',[0.5 0.5 0.5],'LineWidth',2)
    box('off')
    xlim ([0 1000])
    ylim ([-20 20])
    xticks([0:250:1000])
    xticklabels ({'-500','-250','0','250','500'})
    fig.Units = 'centimeters';
    fig.OuterPosition= [10, 10, 10, 10];
    fig.Color='w';
    set(gca,'FontSize',12)
    ylabel('change in beta amplitude (%)');
    xlabel('Time(msec)');
    
    hold on
end

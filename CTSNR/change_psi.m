function [data,baseline]=change_psi(dat1,dat2,dat3,co1,co2,co3);
m=0; plt=[];
for pr=1:size(dat1,1)
    raw_ref1=dat1{pr,1};
    raw_ref2=dat2{pr,1};
    raw_ref3=dat3{pr,1};
    r=0;
    for ct1=1:size(raw_ref1,1)
        for ct2=1:size(raw_ref2,1)
            for ct3=1:size(raw_ref3,1)
                [Pxx_ind,F_i]=mscohere(raw_ref2(ct2,:),raw_ref3(ct3,:),1000,[],1000,1000);
                frange=find(F_i==15):find(F_i==35);
%                  if (sum(Pxx_ind(frange))/sum(Pxx_ind(1:end)))>0.1
                    r=r+1;
                    m=m+1;
                    Pxx_ind_beta=Pxx_ind(frange);
                    filtrange=14+find(Pxx_ind_beta==max(Pxx_ind_beta));
                    [b,a]=butter(2,[(filtrange-5)/(0.5*1000) (filtrange+5)/(0.5*1000)],'bandpass');
                    filt_ref1=filtfilt(b,a,raw_ref1(ct1,:));
                    filt_ref2=filtfilt(b,a,raw_ref2(ct2,:));
                    filt_ref3=filtfilt(b,a,raw_ref3(ct3,:));
                    env_ref1=abs(hilbert(filt_ref1)); env_ref1=smoothdata(env_ref1,'movmean',50);
                    pha_ref2=angle(hilbert(filt_ref2));
                    pha_ref3=angle(hilbert(filt_ref3));
                    
                    dif_pha=wrapToPi(pha_ref2-pha_ref3);
                    
                    [onset1,offset1]=bursts(env_ref1);
                    onset1=horzcat(onset1{:});
                    offset1=horzcat(offset1{:});
                    
                    for o=1:size(onset1,1)
                        %                             onset=onset1{o,1};
                        onset=onset1;
                        el=500;
                        for jj=1:length(onset)
                            if onset(jj)>el && onset(jj)+el<length(dif_pha)
                                epochs_idx(jj,:)=onset(jj)-el:onset(jj)+el;
                                epochs_t(jj,:)=dif_pha(onset(jj)-el:onset(jj)+el);
                            else
                                epochs_idx(jj,1:1001)=NaN;
                                epochs_t(jj,1:1001)=NaN;
                            end
                        end
                        
                        epochs_idx = epochs_idx(any(epochs_idx,2),:);
                        epochs_t = epochs_t(any(epochs_t,2),:);
                        
                        ol=50;
                        for z= 1:size(epochs_idx,1)
                            for w=1:size(epochs_idx,2)
                                if  epochs_idx(z,w)+ol<length(dif_pha)
                                    ep_b1(z,w)=circ_r((dif_pha(epochs_idx(z,w):epochs_idx(z,w)+ol))');
                                    ep_t1(r,w)=circ_r(epochs_t(:,w));
                                end
                            end
                        end
                        ep_b(pr,o,r,:)=zscore(nanmean(ep_b1,1));
                        ep_t(pr,o,r,:)=zscore(nanmean(ep_t1,1));
                        base_t(m,:)=mean(ep_t1(1:el)); 
                        base_b(m,:)=mean(ep_b1(1:el)); 
                        clear ep_b1 ep_t1 onset epochs_idx epochs_t
                    end
                    
                    clear epochs_idx_sur epochs_t_sur
                    for n=1:(length(dif_pha)/1000)
                        idx_sur=randi([el+1,(length(dif_pha)-el)],1,1);
                        epochs_idx_sur(n,:)= idx_sur-el:idx_sur+el;
                        epochs_t_sur(n,:)= dif_pha(idx_sur-el:idx_sur+el);
                    end
                    
                    for z= 1:size(epochs_idx_sur,1)
                        for w=1:size(epochs_idx_sur,2)
                            if  epochs_idx_sur(z,w)+ol<length(dif_pha)
                                ep_b1_s(z,w)=circ_r((dif_pha(epochs_idx_sur(z,w):epochs_idx_sur(z,w)+ol))');
                                ep_t1_s(1,w)=circ_r(epochs_t_sur(:,w));
                            end
                        end
                    end
                    ep_b_s(m,:)=zscore(nanmean(ep_b1_s,1)); clear ep_b1_s
                    ep_t_s(m,:)=zscore(ep_t1_s,1); clear ep_t1_s
%                  end
                clear filt_ref1 filt_ref1 env_ref1  pha_ref2 pha_ref2 onset1 surr_seg dif_pha
                plt=1;
            end
        end
    end
    clear raw_ref1 raw_ref2 raw_ref3
end


if ~isempty(plt)
    cond={'ep_b' 'ep_t'};
    
    for c=1:size(cond,1)
        for cc=1:size(cond,2)
            clear dat1
            dat1=eval(cond{c,cc});
            for bst=1
                %                 bst=1:2
                dum=squeeze(dat1(:,bst,:,300:800)); %%% CHOOSING [-200 300] TO PLOT
                join=[];
                for i=1:size(dum,1)
                    dum1=squeeze(dum(i,:,:));
                    dum2=dum1(any(dum1,2),:);
                    join=[join ; dum2]; clear dum1 dum2
                end
                data{c,cc}(bst,:,:)=join; clear join
            end
        end
    end
    
    data{1,3}=ep_b_s(:,300:800);
    data{1,4}=ep_t_s(:,300:800);
    
    baseline=round([mean(base_b) mean(base_t)],2);
end

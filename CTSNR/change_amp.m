function [data]=change_amp(dat1,dat2)
ref1_amp=[];
ref2_amp=[];
m=0;
for pr=1:size(dat1,1)
    raw_ref1=dat1{pr,1};
    raw_ref2=dat2{pr,1};
    r=0;
    for ct1=1:size(raw_ref1,1)
        for ct2=1:size(raw_ref2,1)
            [Pxx_ind,F_i]=mscohere(raw_ref1(ct1,:),raw_ref2(ct2,:),1000,[],1000,1000);
            frange=find(F_i==15):find(F_i==35);
            
            %             if (sum(Pxx_ind(frange))/sum(Pxx_ind(1:end)))>0.1
            r=r+1;
            m=m+1;
            Pxx_ind_beta=Pxx_ind(frange);
            filtrange=14+find(Pxx_ind_beta==max(Pxx_ind_beta));
            [b,a]=butter(2,[(filtrange-5)/(0.5*1000) (filtrange+5)/(0.5*1000)],'bandpass');
            filt_ref1=filtfilt(b,a,raw_ref1(ct1,:));
            %                  new_d(pr,ct1,:)=filt_ref1;
            % band=[35 60];
            % band=[0.5 3];
            % band=[4,13];
            band=[filtrange-5,filtrange+5];
            
            [b,a]=butter(2,[band(1)/(0.5*1000) (band(end))/(0.5*1000)],'bandpass');
            
            filt_ref2=filtfilt(b,a,raw_ref2(ct2,:));
            %                  new_d(pr,ct2+1,:)=filt_ref2;
            env_ref1=abs(hilbert(filt_ref1)); env_ref1=smoothdata(env_ref1,'movmean',50);
            env_ref2=abs(hilbert(filt_ref2)); env_ref2=smoothdata(env_ref2,'movmean',50);
            
            [onset1,offset1]=bursts(env_ref1);
            onset1=horzcat(onset1{:});
            offset1=horzcat(offset1{:});
            
            for o=1:size(onset1,1)
                %                     onset=onset1{o,1};
                onset=offset1;
                el=500;
                for jj=1:length(onset)
                    if onset(jj)>el && onset(jj)+el<length(env_ref1)
                        ref1_seg(jj,:)=((env_ref1(onset(jj)-el:onset(jj)+el)-median(env_ref1(onset(jj)-el:onset(jj))))./median(env_ref1(onset(jj)-el:onset(jj))));
                        ref2_seg(jj,:)=((env_ref2(onset(jj)-el:onset(jj)+el)-median(env_ref2(onset(jj)-el:onset(jj))))./median(env_ref2(onset(jj)-el:onset(jj))));
                    else
                        ref1_seg(jj,1:1001)=NaN;
                        ref2_seg(jj,1:1001)=NaN;
                    end
                end
                ref1_amp(pr,o,r,:)=nanmedian(ref1_seg,1);
                ref2_amp(pr,o,r,:)=nanmedian(ref2_seg,1);
                clear ref1_seg ref2_seg onset
            end
            for j=1:length(env_ref2)/1000
                idx_sur=randi([501,(length(env_ref2)-el)],1,1);
                surr_seg(j,:)= ((env_ref2(idx_sur-el:idx_sur+el)-median(env_ref2(idx_sur-el:idx_sur)))./median(env_ref2(idx_sur-el:idx_sur)));
            end
            surr_amp(m,:)=nanmedian(surr_seg,1);
        end
        clear filt_ref1 filt_ref1 env_ref1  env_ref2 onset1 surr_seg
    end
end
clear raw_ref1 raw_ref2



cond={'ref1_amp' 'ref2_amp'};

for c=1:size(cond,1)
    for cc=1:size(cond,2)
        clear dat1
        dat1=eval(cond{c,cc});
        for  bst=1
            %             bst=1:2
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

data{1,3}=surr_amp(:,300:800);

end

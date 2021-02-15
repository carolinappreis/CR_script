function [fig]=pha_slip_sub_subref(coh_filts,name)

r=0;
for pr=1:size(coh_filts,1)
    thal=coh_filts{pr,1}(2:end,:);
    
    for ct=1:size(thal,1)
        env=abs(hilbert(thal(ct,:)));
        
        [onset,offset1]=bursts_aligned(env,thal(ct,:));
        if numel(offset1{2,1})>numel(onset{2,1})
            offset1{2,1}=offset1{2,1}(1:length(onset{2,1}));
        end
        dur_long=offset1{2,1}-onset{2,1};
        [dur,dur_idx]=sort(dur_long,'descend');
        on=flip(dur_idx);
        
        r=r+1;
        pha_thal=angle(hilbert(thal(ct,:)));
        non_norm=unwrap(pha_thal);
        non_norm1=diff(non_norm);
        znon_norm=zscore(non_norm1);
        el=400;
        
        for ii=1:length(on)
            if onset{2,1}(on(ii))>el && onset{2,1}(on(ii))+el<length(znon_norm)
                epochs_zi(ii,:)=znon_norm(onset{2,1}(on(ii))-el:onset{2,1}(on(ii))+el);
            end
        end
%         
%         for ii=1:(length(znon_norm)/1000)
%             idx_sur=randi([el+1,(length(znon_norm)-el)],1,1);
%             epochs_zi(ii,:)= znon_norm(idx_sur-el:idx_sur+el);
%         end
%         
%         
        epochs_z=epochs_zi(1:25,:);
        
        for m=1:size(epochs_z,1)
            for ff=1:length(epochs_z(m,:))
                if epochs_z(m,ff)>=1.96 | epochs_z(m,ff)<=-1.96
                    epochs_z1(m,ff)=1;
                else
                    epochs_z1(m,ff)=0;
                end
            end
        end
        epochs_ct(r,:,:)=epochs_z1;
        
        clear pha_thal epochs_z epochs_z1 znon_norm non_norm1 non_norm env onset offset1 dur dur_idx on epochs_zi
    end
    clear thal
end

clearvars -except epochs_ct name coh_filts

load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','squash','blood','sky','aegean');
if name=='bz'
    color_b=blood;
else
    color_b=aegean;
end


sl=squeeze(mean(epochs_ct,1));
bi=1:10:size(sl,2);

for t=1:size(sl,1)
    for r= 1:size(bi,2)
        if r+1<=length(bi)
            ns(t,r)=sum(sl(t,(bi(r):bi(r+1)-1)))./10;
        else
            ns(t,r)=NaN;
        end
    end
end

fig=figure;
subplot(2,1,1)
imagesc(ns)
colorbar
if name=='bz'
    title('BZ contacts')
else
    title('SNR contacts')
end
xlabel ('Time (msec)')
ylabel('Bursts (sorted by length)')
set(gca,'FontSize',12)
xline(40,'r--',{'burst onset'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','LineWidth',2)
xticks([20:20:80])
xlim([20 80])
xticklabels ({'-200','0','200','400'})

pl=squeeze(mean(epochs_ct,2));
pl(pl~=0)=1;
ppl=sum(pl,1)./(size(pl,1));
clear bi;bi=1:10:size(ppl,2);

for r= 1:size(bi,2)
    if r+1<=length(bi)
        ps(1,r)=sum(ppl(1,(bi(r):bi(r+1)-1)))./10;
    else
        ps(1,r)=NaN;
    end
end

subplot(2,1,2)
plot(ps,'Color',color_b,'LineWidth',2)
%%%% ONSET
xline(40,'--',{'burst onset'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','LineWidth',2,'Color',[0.5 0.5 0.5])
xticks([20:20:80])
xlim([20 80])
ylim([0 0.8])
yticks([0:0.2:0.8])
xticklabels ({'-200','0','200','400'})
fig.Units = 'centimeters';
fig.InnerPosition= [10, 10, 14,12];
fig.Color='w';
set(gca,'FontSize',12)
xlabel ('Time (msec)')
ylabel('Probability of phase slip')
box('off')


[stats_pre_pos,p]=ttest(mean(pl(:,200:399),2),mean(pl(:,401:600),2))

end

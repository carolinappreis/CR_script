% function [fig]=pha_slip_ctx_subctx(coh_filts,name)

r=0; tocomp=[]; tocomp_r=[];
for pr=1:size(coh_filts,1)
    ctx=coh_filts{pr,1}(1,:);
    thal=coh_filts{pr,1}(2:end,:);
    
    pha_ctx=angle(hilbert(ctx));
    env=abs(hilbert(ctx));
    
    [onset,offset]=bursts_aligned(env,ctx);
    if numel(offset{2,1})>numel(onset{2,1})
        offset{2,1}=offset{2,1}(1:length(onset{2,1}));
    elseif numel(onset{2,1})>numel(offset{2,1})
        onset{2,1}=onset{2,1}(1:length(offset{2,1}));
    end
    dur_long=offset{2,1}-onset{2,1};
    
    ref=onset{2,1};
    [dur,dur_idx]=sort(dur_long,'descend');
    on=flip(dur_idx(1:10));
    
    for ct=1:size(thal,1)
        r=r+1;
        pha_thal=angle(hilbert(thal(ct,:)));
        
        non_norm=unwrap(pha_ctx)-unwrap(pha_thal); %circdist
        non_norm1=diff(non_norm); 
         non_norm2=(non_norm1-mean(non_norm1))./std(non_norm1);
        znon_norm=zscore(non_norm2);
        el=400;
        start=200; %  duration with no bursts 
        
        for ii=1:length(on)
            if ref(on(ii))>el
                epochs_z(ii,:)=znon_norm(ref(on(ii))-el:ref(on(ii))+el);
                dif_ang(ii,:)=(wrapToPi(non_norm(ref(on(ii))-el:ref(on(ii))+el)));
                after_on(ii,:)=circ_mean((wrapToPi(non_norm(ref(on(ii)):ref(on(ii))+50)))');
            end
        end
        for m=1:size(epochs_z,1)
            if ~isempty(find(abs(epochs_z(m,start:el))>1.96))
                main_slip=(max(find(abs(epochs_z(m,start:el))>1.96)))+start;
                if main_slip+50<length(dif_ang(m,:)) 
                    b_ps(m,1)=circ_mean(wrapToPi((circ_mean((dif_ang(m,main_slip-50:main_slip)'))-after_on(m))));
                    a_ps(m,1)=circ_mean(wrapToPi((circ_mean((dif_ang(m,main_slip:main_slip+50)'))-after_on(m))));
                    bf(m,1)=wrapToPi(circ_mean(dif_ang(m,main_slip-50:main_slip)'));
                    af(m,1)=wrapToPi(circ_mean(dif_ang(m,main_slip:main_slip+50)'));
                else
                    b_ps(m,1)=NaN;
                    a_ps(m,1)=NaN;
                end
            else
                b_ps(m,1)=NaN;
                a_ps(m,1)=NaN;
            end
            
            for ff=1:length(epochs_z(m,:))
                if epochs_z(m,ff)>=1.96 | epochs_z(m,ff)<=-1.96
                    epochs_z1(m,ff)=1;
                else
                    epochs_z1(m,ff)=0;
                end
            end
        end
         if size(epochs_z1,1)==10
        epochs_ct(r,:,:)=epochs_z1;
        tocomp(r,:)=[circ_mean(b_ps(find(~isnan(b_ps)))) circ_mean(a_ps(find(~isnan(a_ps))))];
        tocomp_r(r,:)=[circ_r(b_ps(find(~isnan(b_ps)))) circ_r(a_ps(find(~isnan(a_ps))))];
        abs_r(r,:)=[circ_r(bf(find(~isnan(bf)))) circ_r(af(find(~isnan(af))))];
         end
        clear pha_thal epochs_z epochs_z1 znon_norm non_norm1 non_norm pref_pha b_ps after_on b_ps a_ps
    end
    clear pha_thal onset offset1 dur dur_idx on ctx thal pha_ctx env
end

clearvars -except  name coh_filts epochs_ct tocomp_r tocomp abs_r


load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','squash','blood','sky','aegean');
if name=='bz'
    color_b=blood;
else
    color_b=aegean;
end


fig=figure;
y1=abs_r;
err=std(y1);
x=mean(y1);
bar(x,'EdgeColor','none','FaceColor',color_b,'FaceAlpha',0.4)
hold on
errorbar([1:2],x,err,'.','Color',color_b,'LineWidth',3)
[p,h]=(ttest(y1(:,1),y1(:,2)))
line([1;2],[0.5;0.5],'Color','k');
if h<0.001
    plot([1.4 1.55 1.7] ,[0.55 0.55 0.55],'k*','LineWidth',0.5)
elseif h<0.01
    plot([1.4 1.55] ,[0.55 0.55],'k*','LineWidth',0.5)
elseif h<0.05
    plot(1.5,0.55,'k*','LineWidth',0.5)
end

ylim([0 0.6])
xticklabels({'pre-slip','post-slip'})
ylabel('PSI')
box('off')
set(gca,'FontSize',12)
fig.Units='centimeters';
fig.OuterPosition= [10, 10, 8, 10];
set(fig,'color','w');



maxi=20;

fig=figure;
subplot(1,2,1)
polarhistogram(tocomp(:,1),'FaceColor',color_b,'FaceAlpha',.4,'EdgeColor','none','BinWidth',2*pi/12)
rlim([0 maxi])
set(gca,'FontSize',12)
circ_rtest(tocomp(:,1))


subplot(1,2,2)
polarhistogram(tocomp(:,2),'FaceColor',color_b,'FaceAlpha',.4,'EdgeColor','none','BinWidth',2*pi/12)
rlim([0 maxi])
fig.Units = 'centimeters';
fig.InnerPosition= [10, 10, 14,6];
fig.Color='w';
set(gca,'FontSize',12)
circ_rtest(tocomp(:,2))

fig=figure;
subplot(1,2,1)
% h=histogram(tocomp_r(:,1),'FaceColor',color_b,'FaceAlpha',.4,'EdgeColor','none')
h=histogram(abs_r(:,1),'FaceColor',color_b,'FaceAlpha',.4,'EdgeColor','none')
h.NumBins = 5;
xlim ([0 0.8])
xticks(0:0.2:0.8)
ylim([0 maxi])
yticks([0:15:maxi])
box('off')
xlabel('PSI')
ylabel('Counts')
set(gca,'FontSize',12)

subplot(1,2,2)
% histogram(tocomp_r(:,2),'FaceColor',color_b,'FaceAlpha',.4,'EdgeColor','none')
h=histogram(abs_r(:,2),'FaceColor',color_b,'FaceAlpha',.4,'EdgeColor','none')
h.NumBins = 5;
xlim ([0 0.8])
xticks(0:0.2:0.8)
ylim([0 maxi])
yticks([0:15:maxi])
box('off')
xlabel('PSI')
ylabel('Counts')
fig.Units = 'centimeters';
fig.InnerPosition= [10, 10, 14,6];
fig.Color='w';
set(gca,'FontSize',12)



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
% colorbar
if name=='bz'
    title('CTX-BZ')
else
    title('CTX-SNR')
end
xlabel ('Time (msec)')
ylabel('Bursts (sorted by length)')
set(gca,'FontSize',12)
xline(40,'r--',{'burst onset'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','LineWidth',2)
xticks([20:20:80])
xlim([20 80])
xticklabels ({'-200','0','200','400'})
box('off')

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
% fig=figure;
plot(ps,'Color',color_b,'LineWidth',2)
%%%% ONSET
xline(40,'--',{'burst onset'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','LineWidth',2,'Color','r')
xticks([20:20:80])
xlim([20 80])
ylim([0 0.8])
yticks([0:0.2:0.8])
xticklabels ({'-200','0','200','400'})
fig.Units = 'centimeters';
fig.InnerPosition= [10, 10, 12,12];
fig.Color='w';
set(gca,'FontSize',12)
xlabel ('Time (msec)')
ylabel('Probability of phase slip')
box('off')


[stats_pre_pos,p]=ttest(mean(pl(:,200:399),2),mean(pl(:,401:600),2))

% end

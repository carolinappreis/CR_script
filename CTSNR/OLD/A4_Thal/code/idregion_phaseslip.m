clear all
close all
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal')
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal')

load('BZ_opt.mat');
color_b=[0.5 0 0];
% color_b=[0 0 0.8];
CTX=0;


% for ii=1:length(BZ.idrat)
% data{ii,1}=BZ.filt_thal{BZ.idrat(ii),1}
% end
% data=vertcat(data{:});
bins=[55:50:300];



for ik=1:size(BZ.env_ctx,1)
    clearvars dur
    ref1=BZ.onset_raw_all{ik,1};
    %     ref1_1=BZ.offset_pa_all{ik,1}; %%% OFFSET
    ref1_1=BZ.onset_pa_all{ik,1}; %%% ONSET
    ref2=BZ.offset_raw_all{ik,1};
    if length(ref1) ~= length(ref1_1)
        ref1=ref1(1:length(ref1_1));
        ref2=ref2(1:length(ref1_1));
    end
    %         [dur,dur_idx]=sort(ref2-ref1,'ascend');
    [dur,dur_idx]=sort(ref2-ref1,'ascend');
    %     dur_all(ik,1)=max(dur);
    if max(dur)>=290
        for t=1:length(bins)-1
            r=1;
            for i=1:length(dur)-1
                if (dur(i))>bins(t) && (dur(i))<=bins(t+1)
                    ind_b1{ik,t}(1,r)=ref1_1(dur_idx(i));
                    ind_d1{ik,t}(1,r)=dur(i);
                    r=r+1;
                end
            end
        end
    end
end

m=1;
new_idx=[];
for i =1:size(ind_b1,1)
    if ~isempty(ind_b1{i,1})
        new_idx=[new_idx i];
        for g=1:4
            ind_b1_1{m,g}=ind_b1{i,g};
            ind_d1_1{m,g}=ind_d1{i,g};
        end
        m=m+1;
    end
    
end
BZ.idrat=new_idx;
r= min(cellfun(@length,ind_b1_1));



epochs_zd=[];
for bi=1:size(ind_b1,2)
    dim=[];
    for ik=1:length(BZ.idrat)
        clear ref3 epochs_ct
        ref3=ind_b1_1{ik,bi}(1,1:r(bi));
        
        if CTX==0
            what=size(BZ.phase_thal{BZ.idrat(ik),1},1);
        else
            what=1;
        end
        
        for ct=1:what
            clear epochs_z epochs_z1
            
            if CTX==0
                non_norm=unwrap(BZ.phase_thal{BZ.idrat(ik),1}(ct,:));%circdist
            else
                non_norm=unwrap(BZ.phase_ctx(BZ.idrat(ik),:));
            end
            
            non_norm1=diff(non_norm);
            znon_norm=zscore(non_norm1);
            el=400;
            for ii=1:length(ref3)
                if ref3(ii)>el
                    epochs_z(ii,:)=znon_norm(ref3(ii)-el:ref3(ii)+el);
                end
            end
            for ii=1:size(epochs_z,1)
                for ff=1:length(epochs_z(ii,:))
                    if epochs_z(ii,ff)>=1.96 | epochs_z(ii,ff)<=-1.96
                        epochs_z1(ii,ff)=1;
                    else
                        epochs_z1(ii,ff)=0;
                    end
                end
            end
            
            epochs_zd=[epochs_zd ; epochs_z1];
            epochs_ct(ct,:,:)=epochs_z1;
        end
        epochs_probe(ik,:,:)=squeeze(mean(epochs_ct,1));
        epochs_cr{bi,1}(ik,:,:)=squeeze(mean(epochs_ct,1));
        dim=[dim ; epochs_ct];
        clear epochs_ct
    end
    dim_all{bi,1}=dim;
    slip_b{bi,:}=squeeze(mean(epochs_probe,1));
    clear epochs_probe
end

%% HC's


clear pl ppl ps sl
for f=1:size(dim,1)
    epoch_subj(f,:,:)=[dim_all{1,1}(f,:,:) dim_all{2,1}(f,:,:) dim_all{3,1}(f,:,:) dim_all{4,1}(f,:,:)];
end

sl=squeeze(mean(epoch_subj,1));
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

fig=figure(2)
subplot(2,1,1)
imagesc(ns)
title('BZ')
xlabel ('Time (msec)')
ylabel('Bursts (sorted by length)')
yticks([1 65])
yticklabels({'55','300'})
set(gca,'FontSize',12)
xline(40,'r--',{'burst onset'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','LineWidth',2)
xticks([20:20:80])
xlim([20 80])
xticklabels ({'-200','0','200','400'})
    
pl=squeeze(mean(epoch_subj,2));
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
xticklabels ({'-200','0','200','400'})
fig.Units = 'centimeters';
fig.InnerPosition= [10, 10, 10,12];
fig.Color='w';
set(gca,'FontSize',12)
xlabel ('Time (msec)')
ylabel('Probability of phase slip')
box('off')


[stats_pre_pos,p]=ttest(mean(pl(:,200:399),2),mean(pl(:,401:600),2))


%% POSTER VERSION  ---- avg across probes mat=bursts*time
clear pl ppl ps

sl=cell2mat(slip_b);

% figure()
% imagesc(sl)
% xlabel ('Time(msec)')
% ylabel('Bursts # (Sorted by length)')
%%%%% ONSET
% xticks([200:200:800])
% xlim([200 800])
% xticklabels ({'-200','0','200','400'})
% title('BZ aligned to burst onset')
%%%% OFFSET
% xticks([0:200:600])
% xlim([0 600])
% xticklabels ({'-400','-200','0','200'})
% title('BZ aligned to burst offset')

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



fig=figure()
subplot(2,1,1)
imagesc(ns)
title('BZ-CTX')
xlabel ('Time (msec)')
ylabel('Bursts (sorted by length)')
yticks([1 65])
yticklabels({'55','300'})
% fig.Units = 'centimeters';
% fig.OuterPosition= [10, 10, 10, 10];
% fig.Color='w';
set(gca,'FontSize',12)

%%%% ONSET
xline(40,'r--',{'burst onset'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','LineWidth',2)
xticks([20:20:80])
xlim([20 80])
xticklabels ({'-200','0','200','400'})

%%%% OFFSET
% xline(40,'r--',{'burst onset'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','LineWidth',2)
% xticks([0:20:60])
% xlim([0 60])
% xticklabels ({'-400','-200','0','200'})


%%%%% probability of phase slip


pl=sl;pl(pl~=0)=1;
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
xticklabels ({'-200','0','200','400'})

%%%% OFFSET
% xline(40,'--',{'burst onset'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','LineWidth',2,'Color',[0.5 0.5 0.5])
% xticks([0:20:60])
% xlim([0 60])
% xticklabels ({'-400','-200','0','200'})


fig.Units = 'centimeters';
fig.InnerPosition= [10, 10, 10,12];
fig.Color='w';
set(gca,'FontSize',12)
xlabel ('Time (msec)')
ylabel('Probability of phase slip')
box('off')

[stats_pre_pos,p]=ttest(mean(sl(:,200:399),2),mean(sl(:,401:600),2))


%% CR's NEW 

clear pl ppl ps
for f=1:length(new_idx)
    epoch_subj(f,:,:)=[epochs_cr{1,1}(f,:,:) epochs_cr{2,1}(f,:,:) epochs_cr{3,1}(f,:,:) epochs_cr{4,1}(f,:,:)];
end
% % % (numel((squeeze(mean(epoch_subj,1))==sl)==1))==(size(epoch_subj,2)*size(epoch_subj,3))

sl=squeeze(mean(epoch_subj,1));
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

fig=figure(2)
subplot(2,1,1)
imagesc(ns)
title('BZ-CTX')
xlabel ('Time (msec)')
ylabel('Bursts (sorted by length)')
yticks([1 65])
yticklabels({'55','300'})
set(gca,'FontSize',12)
xline(40,'r--',{'burst onset'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','LineWidth',2)
xticks([20:20:80])
xlim([20 80])
xticklabels ({'-200','0','200','400'})
    
pl=sl;pl(pl~=0)=1;
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
xticklabels ({'-200','0','200','400'})
fig.Units = 'centimeters';
fig.InnerPosition= [10, 10, 10,12];
fig.Color='w';
set(gca,'FontSize',12)
xlabel ('Time (msec)')
ylabel('Probability of phase slip')
box('off')


[stats_pre_pos,p]=ttest(mean(sl(:,200:399),2),mean(sl(:,401:600),2))



%% id
clear pl ppl ps


for f=1:length(new_idx)
    epoch_subj(f,:,:)=[epochs_cr{1,1}(f,:,:) epochs_cr{2,1}(f,:,:) epochs_cr{3,1}(f,:,:) epochs_cr{4,1}(f,:,:)];
end

for iii=1:length(new_idx)
    pl=squeeze(epoch_subj(iii,:,:));
    pl(pl~=0)=1;
    ppl=sum(pl,1)./(size(pl,1));
    clear bi;bi=1:10:size(ppl,2);
    for r= 1:size(bi,2)
        if r+1<=length(bi)
            ps(iii,r)=sum(ppl(1,(bi(r):bi(r+1)-1)))./10;
        else
            ps(iii,r)=NaN;
        end
    end
    clear pl
end
fig=figure
plot(mean(ps),'Color',color_b,'LineWidth',2)
hold on
xline(40,'--',{'burst onset'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','LineWidth',2,'Color',[0.5 0.5 0.5])
xticks([20:20:80])
xlim([20 80])
xticklabels ({'-200','0','200','400'})
fig.Units = 'centimeters';
fig.InnerPosition= [10, 10, 10,12];
fig.Color='w';
set(gca,'FontSize',12)
xlabel ('Time (msec)')
ylabel('Probability of phase slip')
box('off')

[stats_pre_pos,p]=ttest(mean(ps(:,20:39),2),mean(ps(:,41:60),2))

dum=squeeze(mean(epoch_subj,1));
[stats_pre_pos1,p1]=ttest(mean(dum(:,200:399),2),mean(dum(:,401:600),2))
% 
% 
% 
% %%% individual heat plot
% for hh=1:size(epoch_subj,1)
%     clear sl
%     sl=squeeze(epoch_subj(hh,:,:));
%     bi=1:10:size(sl,2);
%     
%     for t=1:size(sl,1)
%         for r= 1:size(bi,2)
%             if r+1<=length(bi)
%                 ns(t,r)=sum(sl(t,(bi(r):bi(r+1)-1)));
%             else
%                 ns(t,r)=NaN;
%             end
%         end
%     end
%     
%     subplot(2,5,hh)
%     imagesc(ns)
%     title('BZ-CTX')
%     xlabel ('Time (msec)')
%     ylabel('Bursts (sorted by length)')
%     yticks([1 65])
%     yticklabels({'55','300'})
%     % fig.Units = 'centimeters';
%     % fig.OuterPosition= [10, 10, 10, 10];
%     % fig.Color='w';
%     set(gca,'FontSize',12)
%     
%     %%%% ONSET
%     xline(40,'r--',{'burst onset'},'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','LineWidth',2)
%     xticks([20:20:80])
%     xlim([20 80])
%     xticklabels ({'-200','0','200','400'})
%     
% end

%a)find the contacts of each probe that are not beta coherent with EEG and
%exclude them
%b)beta filter the coherent contacts 
%c)z-score the beta-filtered 
%d)take the envelope the beta-filtered-z-score
%e)select the contacts that in a window of 200msec before and after the
%phase-aligned-burst-onset have a max(env_z-score)>1,96 
%f) take the median and 25/75prctile of max(env_zcore) across contacts,
%env_zcore across contacts and z_score_beta_filtred 

clear all

cd ('\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A2_Thal\mat')
load 'lesion_param_2000_SNr_new.mat'

time_seg=-2000:1:2000;
time_seg2=-200:1:200;
% figure(1)
% for j=1:size(evoked_all,1)
%     for i=1:size(evoked_all{j},1)
%         subplot(size(evoked_all{j},1),1,i)
%         plot(time_seg,evoked_all{j}(i,:))
%         hold on
%     end
%     close all
% end

beta_coherent=beta_coh_all;

for j=1:size(beta_coh_all,1)
    if ~isempty (beta_coh_all{j});
        for k=1:length(beta_coh_all{j,1})
            if beta_coherent{j}(1,k)<0.1;
                beta_coherent{j}(1,k)=NaN;
            end
        end
    end
end

for j=1:size(beta_coh_all,1)
    if ~isempty (beta_coh_all{j});
        index_beta{j,:}=find (~isnan(beta_coherent{j}(1,:)));
    end
end

figure(1)
for j=1:size(evoked_all,1)
    for i=1:size(index_beta{j},2)
        plot(time_seg,evoked_all{j}(i,:))
        hold on
    end
    close all
end

coh_evoked=[];
for j=1:size(evoked_all,1)
    for i=1:size(index_beta{j},2)
        p=index_beta{j}(i);
        coh_evoked{j,1}(i,:)=evoked_all{j}(p,:);
    end
end
coh_evoked1=coh_evoked(~cellfun('isempty',coh_evoked));
coh_evoked=cell2mat(coh_evoked1);
% plot(mean(coh_evoked));

sum_cvoked=[];
for j=1:size(coh_evoked1,1)
    if size(coh_evoked1{j},1)>1
    sum_cvoked{j,1}(:,1:4001)=sum(coh_evoked1{j});
    else
     sum_cvoked{j,1}(:,1:4001)=(coh_evoked1{j});    
    end
end
sum_cvoked=cell2mat(sum_cvoked);
% plot(sum_cvoked');


samprate=1000;
filt_evoked=[];
for j=1:size(evoked_all,1)
    for i=1:size(index_beta{j},2)
        n=index_beta{j}(i)
        filtrange_1=find(power_all{j}(n,:)== max(power_all{j}(n,15:36)));
        [b,a]=butter(2,[(filtrange_1-5)/(0.5*samprate) (filtrange_1+5)/(0.5*samprate)],'bandpass');
        filt_evoked{j}(i,:)=filtfilt(b,a,evoked_all{j}(n,:));
    end
end


filt_evoked=filt_evoked';


for j=1:size(filt_evoked,1)
    for i=1:size(filt_evoked{j},1)
        zscore_thal{j}(i,:)=zscore(filt_evoked{j}(i,:));
    end
end
zscore_thal=zscore_thal';

for j=1:size(zscore_thal,1)
    for i=1:size(filt_evoked{j},1)
        env_zscore_thal{j}(i,:)=abs(hilbert(zscore_thal{j}(i,:)));
    end
end
env_zscore_thal=env_zscore_thal';

select=[round((length(time_seg))./2-200):round((length(time_seg))./2+200)];

for j=1:size(env_zscore_thal,1)
    for i=1:size(filt_evoked{j},1)
        zselect{j}(i,:)=env_zscore_thal{j}(i,select);
    end
end
zselect=zselect';

for j=1:size(env_zscore_thal,1)
    if ~isempty (zselect{j});
        for i=1:size(zselect{j},1)
            if max(zselect{j}(i,:))<=1.96
                zselect{j}(i,:)=NaN;
            end
        end
    end
end


for j=1:size(zselect,1)
    if ~isempty (zselect{j});
        zbest_idx{j,:}=find (~isnan(zselect{j}(:,i)));
    end
end

zbest=cell(1,1);
filt_best=cell(1,1);
for j=1:size(zselect,1)
    if ~isempty (zbest_idx{j});
        for i =1:size(zbest_idx{j},1);
            n=zbest_idx{j}(i);
            zbest{j,:}(i,:)=zselect{j}(n,:);
            filt_best{j,:}(i,:)=filt_evoked{j}(n,:);
            zscore_best{j,:}(i,:)=zscore_thal{j}(n,:);
        end
    end
end

zbest=cell2mat(zbest);
filt_best=cell2mat(filt_best);
zscore_best=cell2mat(zscore_best);
zbest_stats=prctile(zbest,[25 50 75]);
filt_best_stats=prctile(filt_best,[25 50 75]);
zscore_best_stats=prctile(zscore_best,[25 50 75]);

peak=zeros(size(zbest,1),1);
for i = 1:size(zbest,1)
    peak(i,1)= find(zbest(i,:)==max(zbest(i,:)));
end
stats_peak= prctile(peak,[25 50 75]);

% cd('\Users\wolf5173\Documents\Carolina_code\codes_thal\A2_Thal\mat')
% clearvars -except peak select zbest zbest_stats filt_best filt_best_stats stats_peak zscore_best zscore_best_stats
% save 'A2_2_BZ_surr'





% time_seg2=-200:1:200;
% figure(1)
% boxplot(stats_peak)
%
% figure(2)
% for i = 1:size(zbest,1)
% plot(time_seg2,zbest(i,:))
% hold on
% end
%
% figure(3)
% for i=1:3
%     plot(time_seg2,zbest_stats(i,:))
%     hold on
% end
%
% figure(4)
% for i=1:3
%     plot(time_seg2,filt_best_stats(i,select))
%     hold on
% end
%
% figure(5)
% for i=1:size(zbest,1)
% plot(time_seg2(peak(i)),zbest(peak(i)),'ro')
% hold on
% end
% xlim([-200 200])


%%CHECKS
%
% figure(7)
% for i = 1:size(filt_best,1)
%     subplot(2,1,1)
%     plot(time_seg2,filt_best(i,select))
%     hold on
%     subplot(2,1,2)
%     plot(time_seg2,zbest(i,:))
%     hold on
% 
% end
%

% figure(8)
% for i = 1:size(filt_best,1)
%     subplot(2,1,1)
%     dummy=filt_best(i,:);
%     test1=zscore(abs(hilbert(dummy)));
%     test2=max(test1(1,select));
%     star=find(test1==test2);
%     plot(time_seg(star),dummy(star),'r*');
%     hold on
%     plot(time_seg2,dummy(select));
%     subplot(2,1,2)
%     plot(time_seg2,zbest(i,:))
%     hold on
%     star1=find(zbest(i,:)==max(zbest(i,:)));
%     plot(time_seg2(star1),(zbest(i,(star1))),'r*');
%     time_seg(star)==time_seg2(star1)
%     close all
% end
% %
%
% figure(9)
% subplot(2,1,1)
% r=median(filt_best);
% plot(time_seg2,r(1,select))
% hold on
% test1=zscore(abs(hilbert(r)));
% test2=max(test1(1,select));
% star=find(test1==test2);
% plot(time_seg(star),r(star),'r*');
% subplot(2,1,2)
% s=median(zbest);
% plot(time_seg2,s)
% hold on
% star1=find(s==max(s));
% plot(time_seg2(star1),(s(1,(star1))),'r*');
% time_seg(star)==time_seg2(star1)
% 
% figure(10)
% subplot(2,1,1)
% r=median(filt_best);
% plot(time_seg2,r(1,select))
% hold on
% test1=zscore(abs(hilbert(r)));
% test2=max(test1(1,select));
% star2=find(test1==test2);
% plot(time_seg(star2),r(star2),'r*');
% subplot(2,1,2)
% plot(time_seg2,test1(1,select))
% hold on
% plot(time_seg(star2),test1((star2)),'r*')


% for i =1:size(coh_evoked1{1,1},1)
% subplot(size(coh_evoked1{1,1},1),1,i)
% plot(time_seg,coh_evoked1{1,1}(i,:))
% yticks ([]);
% end
% 
% for i=1:size(filt_evoked{3,1},1)
%     subplot(size(filt_evoked{3,1},1),1,i)
% timea=[-2000:2000];
% plot(timea,filt_evoked{3,1}(i,:))
% yticks ([]);
% end
% 
% 
% for i=1:size(zbest{3,1},1)
%     subplot(size(zbest{3,1},1),1,i)
% timeb=[-200:200];
% plot(timeb,zbest{3,1}(i,:))
% yticks ([]);
% end
% 
% for i=1:size(filt_evoked{3,1},1)
% subplot(size(filt_evoked{3,1},1),1,i)
% timea=[-2000:2000];
% timeb=[-200:200];
% plot(timea,zscore_best{3,1}(i,:))
% hold on
% plot(timeb,zbest{3,1}(i,:), 'LineWidth', 1.5)
% ylim ([-7 7])
% end
% title ('Zscore of beta evoked', 'FontSize',16)
% xlabel ('Time (msec)', 'FontSize', 14)
% ylabel ('Amplitude (uV)', 'Fontsize', 14)


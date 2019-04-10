clear all
cd('C:\Users\wolf5173\Documents\Carolina_code\codes_thal\A2_Thal\mat')

load 'lesion_param_2000_BZ.mat'
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

% figure(1)
% for j=1:size(evoked_all,1)
%     for i=1:size(index_beta{j},2)
%         plot(time_seg,evoked_all{j}(i,:))
%         hold on
%     end
%     close all
% end

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

%
filt_evoked=filt_evoked';

% figure(2)
% for j=1:size(filt_evoked,1)
%     if size(filt_evoked{j},1)~=0
%         for i=1:size(filt_evoked{j},1)
%             subplot(size(filt_evoked{j},1),1,i)
%             plot(time_seg,filt_evoked{j}(i,:))
%         end
%     end
% end


for j=1:size(filt_evoked,1)
    for i=1:size(filt_evoked{j},1)
        env_thal{j}(i,:)=abs(hilbert(filt_evoked{j}(i,:)));
    end
end
env_thal=env_thal';


for j=1:size(env_thal,1)
    for i=1:size(filt_evoked{j},1)
        zscore_env{j}(i,:)=zscore(env_thal{j}(i,:));
    end
end
zscore_env=zscore_env';

select=[round((length(time_seg))./2-200):round((length(time_seg))./2+200)];

for j=1:size(env_thal,1)
    for i=1:size(filt_evoked{j},1)
        zselect{j}(i,:)=zscore_env{j}(i,select);
    end
end

zselect=zselect';

for j=1:size(env_thal,1)
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

for j=1:size(zselect,1)
    if ~isempty (zbest_idx{j});
        for i =1:size(zbest_idx{j},1);
            n=zbest_idx{j}(i);
            zbest{j,:}(i,:)=zselect{j}(n,:)
            filt_best{j,:}(i,:)=filt_evoked{j}(n,:)
        end
    end
end

zbest=cell2mat(zbest);
filt_best=cell2mat(filt_best);
zbest_stats=prctile(zbest,[25 50 75]);
filt_best_stats=prctile(filt_best,[25 50 75]);

peak=zeros(size(zbest,1),1);
for i = 1:size(zbest,1)
    peak(i,1)= find(zbest(i,:)==max(zbest(i,:)));
end
stats_peak= prctile(peak,[25 50 75]);

cd('\Users\wolf5173\Documents\Carolina_code\codes_thal\A2_Thal\mat')
clearvars -except peak select zbest zbest_stats filt_best filt_best_stats stats_peak
save 'A2_SNR_surr'

% time_seg2=-200:1:200;
% figure(1)
% boxplot(stats_peak)
% 
% % figure(2)
% % for i = 1:size(zbest,1)
% % plot(time_seg2,zbest(i,:))
% % hold on
% % end
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

figure(8)
for i = 1:size(filt_best,1)
    subplot(2,1,1)
    dummy=filt_best(i,:);
    test1=zscore(abs(hilbert(dummy)));
    test2=max(test1(1,select));
    star=find(test1==test2);
    plot(time_seg(star),dummy(star),'r*');
    hold on
    plot(time_seg2,dummy(select));
    subplot(2,1,2)
    plot(time_seg2,zbest(i,:))
    hold on
    star1=find(zbest(i,:)==max(zbest(i,:)));
    plot(time_seg2(star1),(zbest(i,(star1))),'r*');
    time_seg(star)==time_seg2(star1)
    close all
end
%     

figure(9)
subplot(2,1,1)
r=median(filt_best);
plot(time_seg2,r(1,select))
hold on
test1=zscore(abs(hilbert(r)));
test2=max(test1(1,select));
star=find(test1==test2);
plot(time_seg(star),r(star),'r*');
subplot(2,1,2)
s=median(zbest);
plot(time_seg2,s)
hold on
star1=find(s==max(s));
plot(time_seg2(star1),(s(1,(star1))),'r*');
time_seg(star)==time_seg2(star1)

figure(10)
subplot(2,1,1)
r=median(filt_best);
plot(time_seg2,r(1,select))
hold on
test1=zscore(abs(hilbert(r)));
test2=max(test1(1,select));
star2=find(test1==test2);
plot(time_seg(star2),r(star2),'r*');
subplot(2,1,2)
plot(time_seg2,test1(1,select))
hold on
plot(time_seg(star2),test1((star2)),'r*')




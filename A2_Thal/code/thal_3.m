clear all
cd('/Users/Carolina/Documents/MATLAB/codes_thal')

load 'lesion_param_2000_SNr.mat'
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
for j=1:size(power_all,1)
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

zbest=cell2mat(zselect');

for i=1:size(zbest,1)
    if (max(zbest(i,:)))<=1.96
        zbest(i,:)=NaN;
    end
end


zbest( ~any(zbest,2), : ) = [];

figure(1)
for i = 1:size(zbest,1)
plot(time_seg2,zbest(i,:))
hold on
end


peak=zeros(size(zbest,1),1);
for i = 1:size(zbest,1)
    peak(i,1)= find(zbest(i,:)==max(zbest(i,:)));
end
stats_peak= prctile(peak,[25 50 75]);

figure(2)
boxplot(stats_peak)

figure(3)
plot(time_seg2,prctile(zbest,25))
hold on
plot(time_seg2,prctile(zbest,75))
plot(time_seg2,prctile(zbest,50))

figure(4)
for i=1:size(zbest,1)
plot(time_seg2(peak(i)),zbest(peak(i)),'ro')
hold on
end
xlim([-200 200])



















% for j=1:size(filt_evoked,1)
% for i=1:size(index_beta{j},2)
%    abs_zscore{j}(i,:)=(abs(zscore_env{j}(i,:)));
%   
% end
% end
% abs_zscore=abs_zscore';
% 
% 
% T=1.96;
% for j=1:size(env_thal,1)
%     for i=1:size(index_beta{j},2)
%         subplot(2,1,1)
%         plot(time_seg,abs_zscore{j}(i,:));
%         hold on
%         plot(time_seg,T*ones(length(abs_zscore{j}(i,:)),1));
%         subplot(2,1,2)
%         plot (time_seg,filt_evoked{j}(i,:));
%         hold on
%         sig_idx=find(abs_zscore{j}(i,:)>=T);
%         plot(time_seg(sig_idx),filt_evoked{j}(i,sig_idx),'r.')
%         dummy=7;
%         close all
%     end
% end
% 
% 
% figure(2)
% for j=1:size(filt_evoked,1)
%     if size(index_beta{j},2)~=0
%         for i=1:size(index_beta{j},2)
%             subplot(size(index_beta{j},2),1,i)
%             plot(time_seg,filt_evoked{j}(i,:))
%             hold on
%             sig_idx=find(abs_zscore{j}(i,:)>=T);
%             plot(time_seg(sig_idx),filt_evoked{j}(i,sig_idx),'r.')
%         end
%     else
%         close all
%     end
%     close all
% end

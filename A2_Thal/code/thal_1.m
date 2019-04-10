clear all
cd('/Users/Carolina/Documents/MATLAB/codes_thal')

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

figure(1)
for j=1:size(evoked_all,1)
    for i=1:size(index_beta{j},2)
        plot(time_seg,evoked_all{j}(i,:))
        hold on
    end
    close all
end

samprate=1000;
filt_evoked=[];
for j=1:size(power_all,1)
    for i=1:size(index_beta{j},2)
        filtrange=find(power_all{j}(i,:)== max(power_all{j}(i,15:36)));
        [b,a]=butter(2,[(filtrange-5)/(0.5*samprate) (filtrange+5)/(0.5*samprate)],'bandpass');
        filt_evoked{j}(i,:)=filtfilt(b,a,evoked_all{j}(i,:));
    end
end
%
filt_evoked=filt_evoked';

figure(2)
for j=1:size(filt_evoked,1)
    for i=1:size(index_beta{j},2)
        subplot(size(index_beta{j},2),1,i)
        plot(time_seg,filt_evoked{j}(i,:))
        hold on
    end
    close all
end


for j=1:size(filt_evoked,1)
    for i=1:size(index_beta{j},2)
        env_thal{j}(i,:)=abs(hilbert(filt_evoked{j}(i,:)));
    end
end
env_thal=env_thal';

for j=1:size(env_thal,1)
    for i=1:size(index_beta{j},2)
        zscore_env{j}(i,:)=zscore(env_thal{j}(i,:));
    end
end
zscore_env=zscore_env';

select=[round((length(time_seg))./2-200):round((length(time_seg))./2+200)];

for j=1:size(env_thal,1)
    for i=1:size(index_beta{j},2)
        zselect{j}(i,:)=zscore_env{j}(i,select);
        seg_filt{j}(i,:)= filt_evoked{j}(i,select);
    end
end
zselect=zselect';
seg_filt=seg_filt';

for j=1:size(env_thal,1)
    for i=1:size(zselect{j},1)
        max(zselect{j}(i,:))
        if max(zselect{j}(i,:))<=1.69
           zselect{j}(i,:)=NaN; 
        end
    end
end


for j=1:size(zselect,1)
    if ~isempty (zselect{j});
        zbest_idx{j,:}=find (~isnan(zselect{j}(:,i)));
    end
end

figure(2)
for j=1:size(seg_filt,1)
    if size(seg_filt{j},1)~=0
        for i=1:size(zbest_idx{j},1)
            fake=(seg_filt{j});
            n=zbest_idx{j};
            subplot(length(n),1,i)
            plot(time_seg2,fake(n(i),:))
        end
    end
end

figure(3)
for j=1:size(seg_filt,1)
    if size(seg_filt{j},1)~=0
        for i=1:size(zbest_idx{j},1)
           p=zbest_idx{j}(i);
            subplot(size(zbest_idx{j},1),1,i)
            plot(time_seg2,seg_filt{j}(p,:))
        end
    end
end

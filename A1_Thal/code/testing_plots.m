

%%CHECK1 - plotting all filtered channel per rat
% for i=1:size(muafilt,1)
%     plot(muafilt(i,:))
%     xlim ([0 500])
% end

%CHECK2 - plotting ctx bursts
% for i=1:length(onset_bursts)
% plot(env(onset_bursts(1,i):offset_bursts(1,i)))
% end

%CHECK3 - plotting envelope with start and ending points
% plot(time,env)
% hold on
% plot(time(begin2),env(begin2),'r.')
% plot(time(ending2),env(ending2),'b.')
% plot(time(length(ending2)),tt(length(ending2)))
% n=1


%%CHECK4 - what is happening in all probes during a specific ctx burst (AVGprobe)
% time_plot=[-500:1:500];
% tt2 (size(time_plot,1):size(time_plot,2))=threshold;
% for i=1:size(output,2)
%      close all
%     subplot(2,1,1)
%     % plot(env(onset_bursts(1,i):offset_bursts(1,i)))
%     plot(time_plot,output_ctx(i,:))
%     xlim ([-500 500])
%     hold on
%     plot(time_plot,tt2)
%     subplot(2,1,2)
%     plot(time_plot,median_probes(i,:)) 
% end

%%CHECK5 - plot median ctx_bursts, median thal_bursts
% time_plot=[-500:1:500];
% subplot(2,1,1)
% plot(time_plot,med_ctx_sub)
% subplot(2,1,2)
% plot(time_plot,med_thal_sub)
% 

 
% figure(1)
% subplot (2,1,1)
% time_plot=[-500:1:500];
% plot(time_plot,median(med_ctx))
% hold on
% plot(time_plot,prctile(med_ctx,75));
% plot(time_plot,prctile(med_ctx,25));
% subplot(2,1,2)
% time_plot=[-500:1:500];
% plot(time_plot,median(med_thal))
% hold on
% plot(time_plot,prctile(med_thal,75));
% plot(time_plot,prctile(med_thal,25));





% 
% 
% figure(1)
% subplot (2,1,1)
% time_plot=[-500:1:500];
% plot(time_plot,med_ctx)
% hold on
% subplot(2,1,2)
% time_plot=[-500:1:500];
% plot(time_plot,med_thal)



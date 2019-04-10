% cd('C:\Users\wolf5173\Documents\Carolina_code\codes_thal\A2_Thal\mat')
% load 'snr_median_zscore.mat'
% time_seg=-2000:1:2000;
% time_seg2=-200:1:200;
% select=[round((length(time_seg))./2-200):round((length(time_seg))./2+200)];
% 
% figure(10)
% subplot(2,1,1)
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
% load 'bz_median_zscore.mat'
% subplot(2,1,1)
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

clearvars -except select filt_best_stats filt_best_stats_bz